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

#include "OSspecific.H"
#include <fstream>
#include <string>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::secondCondMMCcurl<CloudType>::readAxisNames
(
    const dictionary& dict,
    const mmcVarSet& Xi
)
{
    // This stage's own coefficients sub-dictionary. Read here rather than
    // through coeffDict() because the axes are needed in the member
    // initialiser list, before SubModelBase has been constructed.
    const word coeffsName(typeName + word("Coeffs"));

    const dictionary& coeffs = dict.subDict(coeffsName);

    if (coeffs.found("referenceAxes"))
    {
        wordList axes(coeffs.lookup("referenceAxes"));

        if (axes.empty())
        {
            FatalErrorInFunction
                << "referenceAxes in " << typeName << "Coeffs is empty; "
                << "at least one pairing axis is required"
                << exit(FatalError);
        }

        return axes;
    }

    // Default: the modified progress variable followed by the mmcVarSet
    // reference variables, in the order they occupy in the particle's XiR()
    const wordList refNames
    (
        mixParticleModel<CloudType>::orderedRefVarNames(Xi)
    );

    wordList axes(refNames.size() + 1);

    axes[0] = "phiModified";

    forAll(refNames, i)
    {
        axes[i + 1] = refNames[i];
    }

    return axes;
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::resolveAxisIndices
(
    const particleType& p
) const
{
    const wordList& axes = this->refAxisNames();

    axisTableIndex_.setSize(axes.size(), -1);

    const auto& table = p.nameVariableLookUpTable();

    DynamicList<word> missing;

    forAll(axes, i)
    {
        axisTableIndex_[i] = table.indexOf(axes[i]);

        if (axisTableIndex_[i] < 0)
            missing.append(axes[i]);
    }

    if (!missing.empty())
    {
        FatalErrorInFunction
            << "Mixing model " << this->modelType() << " is configured to pair"
            << " on " << axes << nl
            << "but the following are not registered on the particle: "
            << missing << nl
            << "Available: " << table.getAllVarNames()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::secondCondMMCcurl<CloudType>::secondCondMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    // Declaring the pairing axes up front lets the base class read this
    // stage's own normalisation factors and size its own reference space,
    // so nothing here depends on the first conditioning's configuration.
    mixParticleModel<CloudType>
    (
        dict, owner, typeName, Xi, readAxisNames(dict, Xi)
    ),

    CL_(this->coeffDict().lookupOrDefault("CL", 0.5)),

    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),

    meanTimeScale_(this->coeffDict().lookup("meanTimeScale")),

    particleFilter_
    (
        this->coeffDict().template lookupOrDefault<word>
        (
            "particleFilter", "secondCondFlag"
        )
    ),

    filterOnFlag_(particleFilter_ == "secondCondFlag"),

    axisTableIndex_(),

    nPairSamples_
    (
        this->coeffDict().template lookupOrDefault<label>("nPairSamples", 0)
    )
{
    if (particleFilter_ != "secondCondFlag" && particleFilter_ != "none")
    {
        FatalErrorInFunction
            << "Unknown particleFilter " << particleFilter_ << nl
            << "Valid options are: secondCondFlag, none"
            << exit(FatalError);
    }

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
    meanTimeScale_(cm.meanTimeScale_),
    particleFilter_(cm.particleFilter_),
    filterOnFlag_(cm.filterOnFlag_),
    axisTableIndex_(cm.axisTableIndex_),
    nPairSamples_(cm.nPairSamples_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::secondCondMMCcurl<CloudType>::includeParticle
(
    const particleType& p
) const
{
    return !filterOnFlag_ || p.secondCondFlag() == 1;
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::referenceVector
(
    const particleType& p,
    scalarField& xi
) const
{
    // Look the axes up by name from the particle rather than by position in
    // XiR(), so that adding or reordering reference variables in
    // mmcVariablesDefinitions cannot silently change what this stage pairs on.
    if (axisTableIndex_.size() != xi.size())
        resolveAxisIndices(p);

    const auto& table = p.nameVariableLookUpTable();

    forAll(xi, i)
    {
        xi[i] = table.get(axisTableIndex_[i]);
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

    // Note: p.dx() and p.dXiR() are deliberately NOT written here. They are
    // the first conditioning's pair-separation diagnostics, and are sampled as
    // such by MixingPopeCloud's Eulerian statistics; writing them from this
    // stage as well would overwrite whole-cloud values with subset ones and
    // silently mix two different quantities in the same output. This stage
    // reports its own separations in aggregate from reportPairing().

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
void Foam::secondCondMMCcurl<CloudType>::reportPairing() const
{
    // Diagnostics for this stage only.
    //
    // Particle counts:    per-process and global participating totals, plus
    //                     the number of pairs / triples produced.
    // Split-axis usage:   how many times each axis was selected as the k-d
    //                     tree split axis. Reveals whether the tree is
    //                     dominated by the progress variable or by the
    //                     shadow-position coordinates.
    // Pair separation:    mean |dXi| per axis over the pairs actually formed -
    //                     the pairing-quality measure for this stage. Recorded
    //                     here rather than on the particles, which carry the
    //                     first conditioning's dx/dXiR.
    //
    // Every reduce() below is reached on every process: buildParticleList()
    // has no early return, so this hook is collective-safe.

    const wordList& axes = this->refAxisNames();
    const label nAxes = axes.size();

    const label nLocal = this->particleList().size();
    label nGlobal = nLocal;
    reduce(nGlobal, sumOp<label>());

    label nPairs   = 0;
    label nTriples = 0;

    scalarField sumSep(nAxes, 0.0);
    scalarField maxSep(nAxes, 0.0);
    label nSep = 0;

    for (const List<label>& pr : this->particlePairs())
    {
        if (pr.size() == 2) ++nPairs;
        else if (pr.size() == 3) ++nTriples;

        // Separation over every adjacent couple in the group
        for (label k = 1; k < pr.size(); ++k)
        {
            const scalarField& a =
                this->eulerianFieldDataList_[pr[k - 1]].XiR();
            const scalarField& b =
                this->eulerianFieldDataList_[pr[k]].XiR();

            if (a.size() < nAxes || b.size() < nAxes) continue;

            for (label i = 0; i < nAxes; ++i)
            {
                const scalar d = mag(a[i] - b[i]);
                sumSep[i] += d;
                maxSep[i] = max(maxSep[i], d);
            }
            ++nSep;
        }
    }

    label nPairsGlobal   = nPairs;
    label nTriplesGlobal = nTriples;
    label nSepGlobal     = nSep;
    reduce(nPairsGlobal,   sumOp<label>());
    reduce(nTriplesGlobal, sumOp<label>());
    reduce(nSepGlobal,     sumOp<label>());

    List<label> hist = this->splitAxisHistogram();
    hist.setSize(nAxes, 0);

    // Reduce component-wise with the scalar/label overloads. Every process
    // reaches this loop with the same nAxes, so the collective count matches.
    for (label i = 0; i < nAxes; ++i)
    {
        reduce(sumSep[i], sumOp<scalar>());
        reduce(maxSep[i], maxOp<scalar>());

        label h = hist[i];
        reduce(h, sumOp<label>());
        hist[i] = h;
    }

    label totalSplits = 0;
    forAll(hist, i) totalSplits += hist[i];

    Info<< "[" << this->modelType() << "] participating particles: "
        << nGlobal << " global;"
        << " pairs = " << nPairsGlobal
        << ", triples = " << nTriplesGlobal << nl;

    Info<< "[" << this->modelType() << "] split-axis histogram (global):"
        << " total splits = " << totalSplits << nl;

    forAll(hist, i)
    {
        const scalar pct =
            (totalSplits > 0)
              ? 100.0*scalar(hist[i])/scalar(totalSplits)
              : 0.0;

        const scalar meanSep =
            (nSepGlobal > 0) ? sumSep[i]/scalar(nSepGlobal) : 0.0;

        Info<< "    axis " << i << " (" << axes[i] << "): "
            << hist[i] << "  (" << pct << " %)"
            << "   mean |d| = " << meanSep
            << ", max |d| = " << maxSep[i] << nl;
    }
    Info<< endl;

    // Also append one tab-separated row per call to
    //   <case>/postProcessing/<model>.log
    // so the diagnostics can be plotted / grepped without trawling through the
    // solver's main log. Only the master process writes, and it writes to the
    // global case path - runTime.path() would put it inside processor0.
    if (Pstream::master())
    {
        const Time& runTime = this->owner().mesh().time();

        const fileName logDir  = runTime.globalPath()/"postProcessing";
        const fileName logPath = logDir/(this->modelType() + ".log");

        mkDir(logDir);

        const bool writeHeader = !Foam::isFile(logPath);

        std::ofstream os(logPath.c_str(), std::ios::out | std::ios::app);

        if (os.is_open())
        {
            if (writeHeader)
            {
                os << "# time\tnParticles\tnPairs\tnTriples";
                forAll(axes, i) os << '\t' << "split_" << axes[i];
                forAll(axes, i) os << '\t' << "meanSep_" << axes[i];
                os << "\ttotalSplits\n";
            }

            os << runTime.value()
               << '\t' << nGlobal
               << '\t' << nPairsGlobal
               << '\t' << nTriplesGlobal;

            forAll(hist, i) os << '\t' << hist[i];

            forAll(axes, i)
            {
                os << '\t'
                   << ((nSepGlobal > 0) ? sumSep[i]/scalar(nSepGlobal) : 0.0);
            }

            os << '\t' << totalSplits << '\n';
        }
    }
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::SmixList()
{
    const Time& runTime = this->owner().mesh().time();

    // Nothing to sample: behave exactly like the base model.
    // Both conditions are uniform across processes, so the collective calls
    // further down are reached by every rank or by none.
    if (nPairSamples_ <= 0 || !runTime.writeTime())
    {
        mixParticleModel<CloudType>::SmixList();
        return;
    }

    // Number of scalars per logged row - keep in step with writePairSamples()
    const label nCols = 12;

    // ---- 1. enumerate the couples that will actually be mixed --------------
    // A group of three is mixed as two overlapping couples, (0,1) then (1,2),
    // so the unit of interest is the couple, not the group. Couples with a
    // remote member are skipped: they are mixed on the other rank too, and the
    // copy held here is discarded, so its post-mix temperature is meaningless.
    DynamicList<label> coupleP;
    DynamicList<label> coupleQ;

    for (const List<label>& pr : this->particlePairs_)
    {
        for (label k = 1; k < pr.size(); ++k)
        {
            const eulerianFieldData& a = this->eulerianFieldDataList_[pr[k-1]];
            const eulerianFieldData& b = this->eulerianFieldDataList_[pr[k]];

            if (a.local() && b.local())
            {
                coupleP.append(pr[k-1]);
                coupleQ.append(pr[k]);
            }
        }
    }

    label nCouplesGlobal = coupleP.size();
    reduce(nCouplesGlobal, sumOp<label>());

    // Uniform sampling probability, so every couple in the run has the same
    // chance of being logged regardless of which rank holds it.
    const scalar pSample =
        (nCouplesGlobal > 0)
      ? min(1.0, scalar(nPairSamples_)/scalar(nCouplesGlobal))
      : 0.0;

    labelList sampleP(coupleP.size());
    labelList sampleQ(coupleQ.size());
    label nSampled = 0;

    forAll(coupleP, i)
    {
        if (this->owner().rndGen().Random() < pSample)
        {
            sampleP[nSampled] = coupleP[i];
            sampleQ[nSampled] = coupleQ[i];
            nSampled++;
        }
    }

    sampleP.setSize(nSampled);
    sampleQ.setSize(nSampled);

    // ---- 2. state before mixing --------------------------------------------
    const wordList& axes = this->refAxisNames();

    label iPhi = -1;
    forAll(axes, i)
    {
        if (axes[i] == "phiModified") iPhi = i;
    }

    List<scalar> rows(nCols*nSampled, 0.0);

    forAll(sampleP, s)
    {
        const eulerianFieldData& ea = this->eulerianFieldDataList_[sampleP[s]];
        const eulerianFieldData& eb = this->eulerianFieldDataList_[sampleQ[s]];

        const particleType& pa = this->particleList_[ea.particleIndex()];
        const particleType& pb = this->particleList_[eb.particleIndex()];

        // Separation on the progress-variable axis, and over the remaining
        // reference axes (the shadow position for the default axis set)
        const scalar phiA = (iPhi >= 0) ? ea.XiR()[iPhi] : 0.0;
        const scalar phiB = (iPhi >= 0) ? eb.XiR()[iPhi] : 0.0;

        scalar dRefSqr = 0.0;
        forAll(axes, i)
        {
            if (i != iPhi) dRefSqr += sqr(ea.XiR()[i] - eb.XiR()[i]);
        }

        const label r = s*nCols;

        rows[r + 0] = scalar(Pstream::myProcNo());
        rows[r + 1] = scalar(pa.secondCondFlag());
        rows[r + 2] = scalar(pb.secondCondFlag());
        rows[r + 3] = phiA;
        rows[r + 4] = phiB;
        rows[r + 5] = mag(phiA - phiB);
        rows[r + 6] = Foam::sqrt(dRefSqr);
        rows[r + 7] = mag(ea.position() - eb.position());
        rows[r + 8] = pa.T();
        rows[r + 9] = pb.T();
        // columns 10, 11 are the post-mix temperatures, filled below
    }

    // ---- 3. the mixing itself ----------------------------------------------
    mixParticleModel<CloudType>::SmixList();

    // ---- 4. state after mixing ---------------------------------------------
    forAll(sampleP, s)
    {
        const eulerianFieldData& ea = this->eulerianFieldDataList_[sampleP[s]];
        const eulerianFieldData& eb = this->eulerianFieldDataList_[sampleQ[s]];

        const particleType& pa = this->particleList_[ea.particleIndex()];
        const particleType& pb = this->particleList_[eb.particleIndex()];

        rows[s*nCols + 10] = pa.T();
        rows[s*nCols + 11] = pb.T();
    }

    writePairSamples(rows, nCols, nCouplesGlobal);
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::writePairSamples
(
    const List<scalar>& rows,
    const label nCols,
    const label nCouplesGlobal
) const
{
    const Time& runTime = this->owner().mesh().time();

    // Collective: every rank contributes its rows, the master writes one file
    List<List<scalar>> allRows(Pstream::nProcs());
    allRows[Pstream::myProcNo()] = rows;
    Pstream::gatherList(allRows);

    if (!Pstream::master()) return;

    const wordList& axes = this->refAxisNames();

    label nRows = 0;
    forAll(allRows, i) nRows += allRows[i].size()/nCols;

    const fileName dir =
        runTime.globalPath()/"postProcessing"/"secondCondPairs";

    mkDir(dir);

    const fileName fName =
        dir/("secondCondPairs_" + runTime.timeName() + ".dat");

    std::ofstream os(fName.c_str());

    if (!os.is_open())
    {
        WarningInFunction
            << "Could not open " << fName << " for writing" << endl;
        return;
    }

    // Which axes went into each distance column
    DynamicList<word> refAxes;
    forAll(axes, i)
    {
        if (axes[i] != "phiModified") refAxes.append(axes[i]);
    }

    // Foam containers have no operator<< for std::ostream - only for
    // Foam::Ostream - so flatten the axis names into a plain string. A single
    // word is fine as-is, since word derives from std::string.
    auto joinNames = [](const UList<word>& names)
    {
        std::string joined;

        forAll(names, i)
        {
            if (i) joined += ' ';
            joined += names[i];
        }

        return joined;
    };

    os << "# " << this->modelType() << " second-conditioning pair samples\n"
       << "# time            " << runTime.timeName() << "\n"
       << "# pairing axes    " << joinNames(axes) << "\n"
       << "# sampled         " << nRows << " of " << nCouplesGlobal
       << " couples globally (target " << nPairSamples_ << ")\n"
       << "#\n"
       << "# dPhiMod   |difference| on the phiModified axis\n"
       << "# dShadow   Euclidean |difference| over " << joinNames(refAxes)
       << "\n"
       << "# dPhys     |difference| in physical space [m]\n"
       << "# T_*_pre   temperature as the couple was paired [K]\n"
       << "# T_*_post  temperature after mixSpeciesOnly [K]\n"
       << "# flag_*    secondCondFlag; always 1 while particleFilter is"
       << " secondCondFlag\n"
       << "#\n"
       << "# proc\tflag_p\tflag_q\tphiMod_p\tphiMod_q\tdPhiMod\tdShadow"
       << "\tdPhys\tT_p_pre\tT_q_pre\tT_p_post\tT_q_post\n";

    os.setf(std::ios::scientific);
    os.precision(6);

    forAll(allRows, procI)
    {
        const List<scalar>& r = allRows[procI];

        for (label i = 0; i + nCols <= r.size(); i += nCols)
        {
            os << label(r[i]) << '\t'
               << label(r[i+1]) << '\t'
               << label(r[i+2]);

            for (label c = 3; c < nCols; ++c) os << '\t' << r[i+c];

            os << '\n';
        }
    }

    Info<< "[" << this->modelType() << "] wrote " << nRows
        << " pair samples to " << fName << endl;
}


template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*patch*/,
    label /*patchFace*/,
    particle& /*p*/
)
{
    // Not applicable: reference variables are initialised by the
    // first-conditioning model.
    return scalarField(0);
}


template<class CloudType>
const Foam::scalarField Foam::secondCondMMCcurl<CloudType>::XiR0
(
    label /*celli*/,
    particle& /*p*/
)
{
    // Not applicable: reference variables are initialised by the
    // first-conditioning model.
    return scalarField(0);
}


template<class CloudType>
void Foam::secondCondMMCcurl<CloudType>::printInfo()
{
    Info<< "Mixing Model: " << this->modelType() << nl
        << token::TAB << "Reference space: " << this->refAxisNames() << nl
        << token::TAB << "Normalisation:   " << this->Xii_ << nl
        << token::TAB << "Particle filter: " << particleFilter_ << nl
        << token::TAB << "Pairing method:  " << this->pairingMethod_ << nl
        << token::TAB << "Timescale:       aISO" << nl
        << token::TAB << "CL:              " << CL_ << nl
        << token::TAB << "CE:              " << CE_ << nl
        << token::TAB << "meanTimeScale:   " << meanTimeScale_ << nl
        << token::TAB << "Mixes:           Y, T, hA  (NOT phi, XiR, XiC, "
        << "secondCondFlag)" << nl
        << token::TAB << "Pair log:        ";

    if (nPairSamples_ > 0)
    {
        Info<< "up to " << nPairSamples_
            << " couples per write time -> postProcessing/secondCondPairs/"
            << endl;
    }
    else
    {
        Info<< "off (set nPairSamples in " << this->modelType()
            << "Coeffs to enable)" << endl;
    }
}


// ************************************************************************* //
