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
 
    CL_(this->coeffDict().lookupOrDefault("CL", 0.5)),
 
    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),
 
    meanTimeScale_(this->coeffDict().lookup("meanTimeScale"))
{
    // Override the base-class Xii_ with the 4D second-conditioning
    // normalisation. The base populates Xii_ from the first-conditioning
    // XiRNames_ (shadow positions only), which leaves phiModified without
    // its own normaliser and excludes xi_z from the k-d tree splitting.
    // Here we read the 4 explicit keys for the
    // (phiModified, xi_x, xi_y, xi_z) reference space.
    const dictionary& Xim = this->coeffDict().subDict("Xim_i");
    this->Xii_.setSize(4);
    this->Xii_[0] = readScalar(Xim.lookup("phiMod_m"));
    this->Xii_[1] = readScalar(Xim.lookup("sPx_m"));
    this->Xii_[2] = readScalar(Xim.lookup("sPy_m"));
    this->Xii_[3] = readScalar(Xim.lookup("sPz_m"));

    printInfo();
}
 
 
template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const secondCondMMCcurl<CloudType>& cm
)
:
    mixParticleModel<CloudType>(cm),
    CL_(cm.CL_),
    CE_(cm.CE_),
    meanTimeScale_(cm.meanTimeScale_)
{
    // Mirror the primary-constructor override of the base-class Xii_ so
    // a cloned model also carries the 4D normalisation. The base copy
    // constructor re-runs getXiNormalisation() which yields only the 3
    // first-conditioning entries.
    const dictionary& Xim = this->coeffDict().subDict("Xim_i");
    this->Xii_.setSize(4);
    this->Xii_[0] = readScalar(Xim.lookup("phiMod_m"));
    this->Xii_[1] = readScalar(Xim.lookup("sPx_m"));
    this->Xii_[2] = readScalar(Xim.lookup("sPy_m"));
    this->Xii_[3] = readScalar(Xim.lookup("sPz_m"));

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
 
        // 4D reference space for the k-d tree:
        //   XiR[0] = phi_modified  (tight normalisation → dominant split axis)
        //   XiR[1] = xi_x          (shadow position x, from particle XiR()[0])
        //   XiR[2] = xi_y          (shadow position y, from particle XiR()[1])
        //   XiR[3] = xi_z          (shadow position z, from particle XiR()[2])
        //
        // The particle's XiR() field holds the first-conditioning reference
        // variables (shadow positions) set by the base solver.  We read them
        // positionally here; their names are not needed.
        eulerianFields.XiR().resize(4);
        eulerianFields.XiR()[0] = iter().phiModified();
        eulerianFields.XiR()[1] = iter().XiR()[0];
        eulerianFields.XiR()[2] = iter().XiR()[1];
        eulerianFields.XiR()[3] = iter().XiR()[2];
 
        // magSqrRefVar is not used in the aISO timescale; set to 4 zeros
        eulerianFields.magSqrRefVar().resize(4, 0.0);
 
        // Interpolate Eulerian transport properties at the particle location
        eulerianFields.Rand()   = 0.0;
        eulerianFields.DEff()   = DEff_intp_.interpolate(pos, cellI, faceI);
        eulerianFields.D()      = D_intp_   .interpolate(pos, cellI, faceI);
        eulerianFields.Dt()     = Dt_intp_  .interpolate(pos, cellI, faceI);
        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos, cellI, faceI);
        eulerianFields.mu()     = mu_intp_  .interpolate(pos, cellI, faceI);
        eulerianFields.vb()     = vb_intp_  .interpolate(pos, cellI, faceI);
 
        this->eulerianFieldDataList_.append(std::move(eulerianFields));
    }
 
    // Pair the flagged particles locally (no parallel exchange for second
    // conditioning — consistent with the local-pairing-only design note)
    this->findPairs(this->eulerianFieldDataList_, this->particlePairs_);

    // -- Diagnostics ---------------------------------------------------------
    // Particle counts: per-process and global flagged-particle totals,
    //                  plus the number of pairs / triples produced.
    // Split-axis usage: how many times each XiR axis was selected as the
    //                   k-d tree split axis during findPairs(). Reveals
    //                   whether the tree is dominated by phi or by the
    //                   shadow-position coordinates.
    {
        const label nLocal = this->particleList_.size();
        label nGlobal = nLocal;
        reduce(nGlobal, sumOp<label>());

        label nPairs   = 0;
        label nTriples = 0;
        for (const List<label>& pr : this->particlePairs_)
        {
            if (pr.size() == 2) ++nPairs;
            else if (pr.size() == 3) ++nTriples;
        }
        label nPairsGlobal   = nPairs;
        label nTriplesGlobal = nTriples;
        reduce(nPairsGlobal,   sumOp<label>());
        reduce(nTriplesGlobal, sumOp<label>());

        // Sum the per-axis split counts across processors so the histogram
        // reflects the global picture.
        List<label> hist = this->splitAxisHistogram();
        forAll(hist, i)
        {
            label v = hist[i];
            reduce(v, sumOp<label>());
            hist[i] = v;
        }

        label totalSplits = 0;
        forAll(hist, i) totalSplits += hist[i];

        Info<< "[secondCondMMCcurl] flagged particles: "
            << nGlobal << " global"
            << " (local rank0 = " << nLocal << ");"
            << " pairs = " << nPairsGlobal
            << ", triples = " << nTriplesGlobal << nl;

        Info<< "[secondCondMMCcurl] split-axis histogram (global):"
            << " total splits = " << totalSplits << nl;

        forAll(hist, i)
        {
            const word axisName =
                (i == 0) ? word("phiModified")
              : (i == 1) ? word("xi_x")
              : (i == 2) ? word("xi_y")
              : (i == 3) ? word("xi_z")
              :            word("XiR[" + Foam::name(i) + "]");

            const scalar pct =
                (totalSplits > 0)
                  ? 100.0*scalar(hist[i])/scalar(totalSplits)
                  : 0.0;

            Info<< "    axis " << i << " (" << axisName << "): "
                << hist[i] << "  ("
                << pct << " %)" << nl;
        }
        Info<< endl;

	// Also append one tab-separated row per call to
        //   <case>/postProcessing/secondCondMMCcurl.log
        // so the diagnostics can be plotted / grepped without trawling
        // through the solver's main log. Only the master process writes.
        if (Pstream::master())
        {
            const Time& runTime = this->owner().mesh().time();

            const fileName logDir  = runTime.path()/"postProcessing";
            const fileName logPath = logDir/"secondCondMMCcurl.log";

            mkDir(logDir);

            // Header on first write only (empty / missing file)
            const bool writeHeader = !Foam::isFile(logPath);

            std::ofstream os
            (
                logPath.c_str(),
                std::ios::out | std::ios::app
            );

            if (os.is_open())
            {
                if (writeHeader)
                {
                    os << "# time\tnFlagged\tnPairs\tnTriples"
                       << "\tphiMod\txi_x\txi_y\txi_z"
                       << "\ttotalSplits\n";
                }

                os << runTime.value()
                   << '\t' << nGlobal
                   << '\t' << nPairsGlobal
                   << '\t' << nTriplesGlobal;

                forAll(hist, i) os << '\t' << hist[i];

                os << '\t' << totalSplits << '\n';
            }
        }
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
 
    if (A > VSMALL)
        tauP = (1.0 / pEulFields.vb())
             * (sqr(pEulFields.DeltaE()) / (CE_ * A));
 
    if (B > VSMALL)
        tauQ = (1.0 / qEulFields.vb())
             * (sqr(qEulFields.DeltaE()) / (CE_ * B));
 
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
 
    const scalar mixExtent = 1.0 - Foam::exp(-deltaT / (tauMix + VSMALL));
 
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
 
    // Mix species-only: Y, T, hA — NOT phi, XiR, XiC, or secondCondFlag
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
 
    // Mix temperature T
    {
        scalar TAv = (p.wt() * p.T() + q.wt() * q.T()) / wtSum;
        p.T() = p.T() + mixExtent * (TAv - p.T());
        q.T() = q.T() + mixExtent * (TAv - q.T());
    }
 
    // Mix species Y
    {
        scalarField YAv = (p.wt() * p.Y() + q.wt() * q.Y()) / wtSum;
        p.Y() = p.Y() + mixExtent * (YAv - p.Y());
        q.Y() = q.Y() + mixExtent * (YAv - q.Y());
    }
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
        << token::TAB << "Reference space: (phiModified, xi_x, xi_y, xi_z)" << nl
        << token::TAB << "Particle filter: secondCondFlag == 1 only" << nl
        << token::TAB << "k-d tree:        splits on XiR (ncond = 3+i, no physical coord)" << nl
        << token::TAB << "Timescale:       aISO" << nl
        << token::TAB << "CL:              " << CL_           << nl
        << token::TAB << "CE:              " << CE_           << nl
        << token::TAB << "meanTimeScale:   " << meanTimeScale_ << nl
        << token::TAB << "Mixes:           Y, T, hA  (NOT phi, XiR, XiC, secondCondFlag)"
        << endl;
}
 
 
// ************************************************************************* //

