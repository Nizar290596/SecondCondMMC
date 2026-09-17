/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
 
    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.
 
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
 
    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 
\*---------------------------------------------------------------------------*/
 
// Note: this file is included by secondCondMMCcurl.H via NoRepository.
// interpolationCellPoint is available through the mixParticleModel.H include chain.

#include "OSspecific.H"
#include <fstream>
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    mixParticleModel<CloudType>(dict, owner, typeName, Xi),
 
    CE_(this->readMixingConstant()),
    includeShadowPositions_(this->coeffDict().lookupOrDefault("includeShadowPositions", true)),
    diagnosticInterval_(this->coeffDict().lookupOrDefault("diagnosticInterval", label(100))),
 
    meanTimeScale_(this->coeffDict().lookup("meanTimeScale"))
{
    configureReferenceSpace();

    printInfo();
}
 
 
template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const secondCondMMCcurl<CloudType>& cm
)
:
    mixParticleModel<CloudType>(cm),
    CE_(cm.CE_),
    includeShadowPositions_(cm.includeShadowPositions_),
    diagnosticInterval_(cm.diagnosticInterval_),
    meanTimeScale_(cm.meanTimeScale_)
{
    configureReferenceSpace();

    printInfo();
}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::buildParticleList()
{
    // Set up interpolators for Eulerian transport fields
    // (identical to mixParticleModel::buildParticleList())
    interpolationCellPoint<scalar> DEff_intp_(this->DEff_);
    interpolationCellPoint<scalar> D_intp_   (this->D_);
    interpolationCellPoint<scalar> Dt_intp_  (this->Dt_);
    interpolationCellPoint<scalar> DeltaE_intp_(this->DeltaE_);
    interpolationCellPoint<scalar> mu_intp_  (this->mu_);
    interpolationCellPoint<scalar> vb_intp_  (this->vb_);
 
    this->particleList_.clear();
    this->eulerianFieldDataList_.clear();
    this->particlePairs_.clear();
 
    label particleInd = 0;
 
    forAllIters(this->owner(), iter)
    {
        // Filter: include only particles flagged for second conditioning
        if (iter().secondCondFlag() != 1)
            continue;
 
        // Append pointer to the particle (no ownership)
        this->particleList_.append(iter.get());
 
        eulerianFieldData eulerianFields;
 
        const vector pos   = iter().position();
        const label  cellI = iter().cell();
        const label  faceI = iter().face();
 
        eulerianFields.particleIndex()  = particleInd++;
        eulerianFields.processorIndex() = Pstream::myProcNo();
        eulerianFields.position()       = pos;
 
        eulerianFields.XiR().resize(this->Xii_.size());
        eulerianFields.XiR()[0] = iter().phiModified();
        if (includeShadowPositions_)
        {
            const auto& index = this->XiR_.rVarInXiR();
            eulerianFields.XiR()[1] = iter().XiR()[index["sPx"]];
            eulerianFields.XiR()[2] = iter().XiR()[index["sPy"]];
            eulerianFields.XiR()[3] = iter().XiR()[index["sPz"]];
        }
        eulerianFields.magSqrRefVar().resize(this->Xii_.size(), 0.0);

        // Interpolate Eulerian transport properties at the particle location
        eulerianFields.Rand()   = this->owner().rndGen().Random();
        eulerianFields.DEff()   = DEff_intp_.interpolate(pos, cellI, faceI);
        eulerianFields.D()      = D_intp_   .interpolate(pos, cellI, faceI);
        eulerianFields.Dt()     = Dt_intp_  .interpolate(pos, cellI, faceI);
        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos, cellI, faceI);
        eulerianFields.mu()     = mu_intp_  .interpolate(pos, cellI, faceI);
        eulerianFields.vb()     = vb_intp_  .interpolate(pos, cellI, faceI);
 
        this->eulerianFieldDataList_.append(std::move(eulerianFields));
    }
 
    const label nLocal = this->particleList_.size();
    const label nAllLocal = this->owner().size();
    scalar flaggedWeight = 0, totalWeight = 0, phiModifiedSecondMoment = 0;
    forAllIters(this->owner(), particle)
    {
        totalWeight += particle().wt();
        if (particle().secondCondFlag() == 1)
        {
            flaggedWeight += particle().wt();
            phiModifiedSecondMoment += particle().wt()*sqr(particle().phiModified());
        }
    }
    if (Pstream::parRun()
        && this->pairingMethod_.method() != particlePairingMethod::localPairing)
        this->correctParticleListParallel();
    else
        this->findPairs(this->eulerianFieldDataList_, this->particlePairs_);

    if (this->owner().mesh().time().timeIndex() % diagnosticInterval_ == 0)
    {
        label nGlobal = nLocal, participating = 0;
        scalar maxDx = 0, maxDphi = 0;
        reduce(nGlobal, sumOp<label>());
        label nAll = nAllLocal;
        reduce(nAll, sumOp<label>());
        reduce(flaggedWeight, sumOp<scalar>());
        reduce(totalWeight, sumOp<scalar>());
        reduce(phiModifiedSecondMoment, sumOp<scalar>());
        for (const auto& group : this->particlePairs_)
        {
            label firstRank = Pstream::nProcs();
            for (label i : group)
                firstRank = min(firstRank, this->eulerianFieldDataList_[i].processorIndex());
            if (firstRank != Pstream::myProcNo()) continue;
            participating += group.size();
            forAll(group, i) for (label j=i+1; j<group.size(); ++j)
            {
                const auto& p = this->eulerianFieldDataList_[group[i]];
                const auto& q = this->eulerianFieldDataList_[group[j]];
                maxDx = max(maxDx, mag(p.position()-q.position()));
                maxDphi = max(maxDphi, mag(p.XiR()[0]-q.XiR()[0]));
            }
        }
        reduce(participating, sumOp<label>());
        reduce(maxDx, maxOp<scalar>());
        reduce(maxDphi, maxOp<scalar>());
        Info<< "secondCondMMCcurl: unique flagged=" << nGlobal
            << "/" << nAll << ", flagged mass fraction=" << flaggedWeight/max(totalWeight, VSMALL)
            << ", paired=" << participating << ", max pair dx=" << maxDx
            << ", max pair deltaPhiModified=" << maxDphi << nl;
        const scalar tau = this->owner().secondCondTauOU();
        if (tau > 0 && flaggedWeight > VSMALL)
            Info<< "OU reference diffusion: weighted mean D_phiModified="
                << sqr(this->owner().secondCondBeta())/tau*phiModifiedSecondMoment/flaggedWeight
                << " [1/s]; compare with Nphi using a consistent conditioning convention" << nl;
    }
}
 
 
template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::mixpair
(
    particleType& p,
    const eulerianFieldData& pEulFields,
    particleType& q,
    const eulerianFieldData& qEulFields,
    scalar& deltaT
)
{
    if (p.wt() + q.wt() <= 0)
        return;
 
    // aISO timescale — identical formulation to MMCcurl
    //   tau = (1/vb) * DeltaE^2 / (CE * (D + Dt))
    scalar tauP = 1e30;
    scalar tauQ = 1e30;
 
    const scalar A = pEulFields.D() + pEulFields.Dt();
    const scalar B = qEulFields.D() + qEulFields.Dt();
 
    if (A > VSMALL && pEulFields.vb() > VSMALL)
        tauP = (1.0 / pEulFields.vb())
             * (sqr(pEulFields.DeltaE()) / (CE_ * A));
 
    if (B > VSMALL && qEulFields.vb() > VSMALL)
        tauQ = (1.0 / qEulFields.vb())
             * (sqr(qEulFields.DeltaE()) / (CE_ * B));
 
    if (this->mixingTimeScale_ == "prescribed")
        tauP = tauQ = this->prescribedTauMix_;

    if (tauP >= 1e30 || tauQ >= 1e30)
        return;
 
    scalar tauMix = 0.0;
 
    if (meanTimeScale_)
        tauMix = 2.0
           / (
                  1.0 / (tauP + VSMALL)
                + 1.0 / (tauQ + VSMALL)
             );
    else
        tauMix = min(tauP, tauQ);
 
    const scalar mixExtent = this->pairMixingExtent(pEulFields, qEulFields, deltaT, tauMix);
 
    // Set diagnostic distance fields on the particles
    scalar dx_pq = Foam::sqrt
    (
        sqr(pEulFields.position().x() - qEulFields.position().x())
      + sqr(pEulFields.position().y() - qEulFields.position().y())
      + sqr(pEulFields.position().z() - qEulFields.position().z())
    );
    p.dx() = dx_pq;
    q.dx() = dx_pq;
 
    forAll(p.dXiR(), i)
    {
        p.dXiR()[i] = mag(p.XiR()[i] - q.XiR()[i]);
        q.dXiR()[i] = p.dXiR()[i];
    }
 
    // Sparse mixing changes thermochemical state, not the dense reference.
    mixSpeciesOnly(p, q, mixExtent);
}
 
 
template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::mixSpeciesOnly
(
    particleType& p,
    particleType& q,
    scalar mixExtent
)
{
    const scalar wtSum = p.wt() + q.wt();
    if (wtSum < VSMALL)
        return;
 
    // Mix enthalpy hA
    {
        scalar hAv = (p.wt() * p.hA() + q.wt() * q.hA()) / wtSum;
        p.hA() = p.hA() + mixExtent * (hAv - p.hA());
        q.hA() = q.hA() + mixExtent * (hAv - q.hA());
    }
 
    // Mix species Y
    {
        scalarField YAv = (p.wt() * p.Y() + q.wt() * q.Y()) / wtSum;
        p.Y() = p.Y() + mixExtent * (YAv - p.Y());
        q.Y() = q.Y() + mixExtent * (YAv - q.Y());
    }
    // Mixture fraction and other transported coupling scalars must follow
    // species mixing; the dense reference phi and shadow coordinates do not.
    forAll(p.XiC(), i)
    {
        const scalar mean = (p.wt()*p.XiC()[i] + q.wt()*q.XiC()[i])/wtSum;
        p.XiC()[i] += mixExtent*(mean-p.XiC()[i]);
        q.XiC()[i] += mixExtent*(mean-q.XiC()[i]);
    }
    // Temperature is derived from the mixed composition and enthalpy.
    p.T() = this->owner().composition().particleMixture(p.Y()).THa(p.hA(), p.pc(), p.T());
    q.T() = this->owner().composition().particleMixture(q.Y()).THa(q.hA(), q.pc(), q.T());
}
 
 
template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*patch*/,
    label /*patchFace*/,
    particle& /*p*/
)
{
    // Not applicable: secondCondMMCcurl does not initialise boundary XiR values.
    return scalarField(0);
}
 
 
template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*celli*/,
    particle& /*p*/
)
{
    // Not applicable: secondCondMMCcurl does not initialise cell XiR values.
    return scalarField(0);
}
 
 
template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::printInfo()
{
    Info<< "Mixing Model: " << this->modelType() << nl
        << token::TAB << "Reference space: "
        << (includeShadowPositions_ ? "(phiModified,xi,x)" : "(phiModified,x)") << nl
        << token::TAB << "Particle filter: secondCondFlag == 1 only" << nl
        << token::TAB << "Physical localization: " << this->physicalLocalization_ << nl
        << token::TAB << "Pairing: " << this->pairingMethod_ << nl
        << token::TAB << "Timescale: " << this->mixingTimeScale_ << nl
        << token::TAB << "CE:              " << CE_           << nl
        << token::TAB << "meanTimeScale:   " << meanTimeScale_ << nl
        << token::TAB << "Mixes: Y, hA, XiC; reconstructs T; preserves dense reference and selection"
        << endl;
}
 
 
// ************************************************************************* //



template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::configureReferenceSpace()
{
    const dictionary& d = this->coeffDict().subDict("Xim_i");
    const bool legacy = d.found("phiModified_m");
    const scalar phiScale = readScalar(d.lookup(d.found("phiMod_m") ? "phiMod_m" : "phiModified_m"));
    if (legacy && d.found("phiMod_m") && phiScale != readScalar(d.lookup("phiModified_m")))
        FatalErrorInFunction << "Conflicting phiMod_m and phiModified_m" << exit(FatalError);
    if (legacy) WarningInFunction << "phiModified_m is accepted as an alias; use phiMod_m." << nl;
    this->Xii_.setSize(includeShadowPositions_ ? 4 : 1);
    this->Xii_[0] = phiScale;
    if (includeShadowPositions_)
    {
        wordList names(3);
        names[0] = "sPx"; names[1] = "sPy"; names[2] = "sPz";
        forAll(names, i)
        {
            if (!this->XiR_.rVarInXiR().found(names[i]))
                FatalErrorInFunction << "Missing reference " << names[i] << exit(FatalError);
            this->Xii_[i+1] = readScalar(d.lookup(names[i]+"_m"));
        }
    }
    for (scalar scale : this->Xii_)
        if (!(scale > 0 && std::isfinite(scale)))
            FatalErrorInFunction << "All second-conditioning normalizations must be positive and finite" << exit(FatalError);
    if (diagnosticInterval_ < 1)
        FatalErrorInFunction << "diagnosticInterval must be >= 1" << exit(FatalError);
    if (!this->physicalLocalization_)
        WarningInFunction << "Physical localization is disabled; this omits the paper's x-space conditioning." << nl;
    if (this->mixingTimeScale_ != "aISO" && this->mixingTimeScale_ != "prescribed")
        FatalErrorInFunction << "Second conditioning supports aISO or prescribed timescales only" << exit(FatalError);
    Info<< "Second reference normalizations: " << this->Xii_ << nl;
}
