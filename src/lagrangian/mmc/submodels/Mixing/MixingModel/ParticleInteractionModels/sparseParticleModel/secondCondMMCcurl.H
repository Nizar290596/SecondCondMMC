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

Class
    Foam::secondCondMMCcurl

Description
    Second-conditioning mixing model for the MMC framework.

    It is an ordinary sparse mixing model that differs from the first
    conditioning in exactly three respects, each expressed through one of the
    hooks mixParticleModel provides:

      includeParticle()   restricts the stage to the flagged subset
      referenceVector()   pairs in this stage's own reference space
      reportPairing()     writes this stage's diagnostics

    Everything else - Eulerian interpolation, the k-d bisection, the parallel
    exchange - is inherited, so this stage supports the same local, subVolumes
    and global pairing methods as any other.

    The stage is configured entirely by its own coefficients sub-dictionary and
    reads nothing belonging to the first conditioning. In particular the
    pairing axes are named explicitly and are looked up by name from the
    particle at run time, so this model has no positional dependence on the
    order of the mmcVarSet reference variables.

    Which particles the stage acts on is a cloud-level property
    (secondCondFlag, assigned from secondConditioning/R) because the reaction,
    coupling and statistics all need the same subset; the model only states
    which predicate it applies.

    Usage in cloudProperties:
    \verbatim
    subModels
    {
        mixingModel             phiMMCcurl;      // mixes phi for all particles
        secondCondMixingModel   secondCondMMCcurl;

        secondCondMMCcurlCoeffs
        {
            r_i             1e-3;

            // The axes this stage pairs on, in order. Any scalar registered on
            // the particle may be used (phi, phiModified, omegaOU, T, the
            // reference variables, the species, ...). Defaults to
            // phiModified followed by the mmcVarSet reference variables.
            referenceAxes   ( phiModified sPx sPy sPz );

            // One normalisation factor per axis above, keyed "<axis>_m".
            // A tight phiModified_m makes the k-d tree split preferentially on
            // the progress variable; looser values on the shadow positions give
            // spatial context at lower priority.
            Xim_i
            {
                phiModified_m   0.01;
                sPx_m           0.05;
                sPy_m           0.05;
                sPz_m           0.05;
            }

            // Which particles take part. secondCondFlag (default) restricts
            // the stage to the subset; none applies it to the whole cloud.
            particleFilter  secondCondFlag;

            pairingMethod   local;

            // aISO timescale coefficients (same meaning as MMCcurl)
            CE              0.1;
            CL              0.5;
            meanTimeScale   true;
        }
    }
    \endverbatim

    Note on dictionary keys: the normalisation entry for the progress variable
    is keyed by its axis name, so it is "phiModified_m" and not "phiMod_m".

    Only Y, T and hA are mixed between the paired particles. phi, XiR, XiC and
    secondCondFlag are deliberately not mixed - they are, or feed, the
    conditioning coordinates.

\*---------------------------------------------------------------------------*/

#ifndef secondCondMMCcurl_H
#define secondCondMMCcurl_H

#include "mixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class secondCondMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class secondCondMMCcurl
:
    public mixParticleModel<CloudType>
{
    // Typedefs (mirror MMCcurl conventions)
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData
        eulerianFieldData;

    // Private data

        //- aISO timescale model coefficients (same as MMCcurl)
        const scalar CL_;
        const scalar CE_;

        //- true  → harmonic mean of the two particle mixing time scales
        //- false → minimum of the two particle mixing time scales
        const Switch meanTimeScale_;

        //- Name of the particle filter this stage applies
        //  "secondCondFlag" (default) or "none"
        const word particleFilter_;

        //- true when particleFilter_ == "secondCondFlag"
        const bool filterOnFlag_;

        //- Index of each pairing axis in the particle's nameVariableTable.
        //  Resolved on first use - the table is populated when a particle is
        //  constructed, so it cannot be resolved in this constructor. The
        //  index is shared by every particle of the cloud.
        mutable labelList axisTableIndex_;


    // Private member functions

        //- Mix only Y, T, and hA between paired particles.
        //  Deliberately does NOT mix phi, XiR, XiC, or secondCondFlag.
        static void mixSpeciesOnly
        (
            particleType& p,
            particleType& q,
            scalar mixExtent
        );

        //- Resolve axisTableIndex_ from a particle, once
        void resolveAxisIndices(const particleType& p) const;

        //- Read the pairing axis names from this model's own coeffDict,
        //  defaulting to phiModified plus the mmcVarSet reference variables
        static wordList readAxisNames
        (
            const dictionary& dict,
            const mmcVarSet& Xi
        );


protected:

    // Stage customisation hooks

        //- Restrict the stage to the flagged subset
        virtual bool includeParticle(const particleType& p) const override;

        //- This stage's reference-space coordinates, looked up by name
        virtual void referenceVector
        (
            const particleType& p,
            scalarField& xi
        ) const override;

        //- Particle counts, split-axis histogram and pair separations
        virtual void reportPairing() const override;


public:

    //- Runtime type information
    TypeName("secondCondMMCcurl");


    // Constructors

        //- Construct from dictionary
        secondCondMMCcurl
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi
        );

        //- Construct copy
        secondCondMMCcurl(const secondCondMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType>> clone() const
        {
            return autoPtr<CloudMixingModel<CloudType>>
            (
                new secondCondMMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~secondCondMMCcurl() = default;


    // Member Functions

        //- Pair-wise mixing: aISO timescale + mixSpeciesOnly
        virtual void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        ) override;

        //- Boundary initialisation. Reference variables are initialised by the
        //  first-conditioning model (MixingPopeCloud calls mixing().XiR0()),
        //  so this stage never supplies them.
        virtual const scalarField XiR0
        (
            label patch,
            label patchFace,
            particle& p
        ) override;

        virtual const scalarField XiR0(label celli, particle& p) override;

    // Info

        void printInfo();
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "secondCondMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
