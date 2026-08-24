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
 
    Operates exclusively on the particle subset flagged with
    secondCondFlag == 1 and pairs particles in the 4-dimensional reference
    space
 
        (phi_modified, xi_x, xi_y, xi_z)
 
    where
      - phi_modified = phiModified = phi*exp(beta*omegaOU), the modified
                       progress variable recomputed each timestep by
                       MixingPopeCloud::updateOUProcess()
      - xi = shadow positions stored in the particle's existing XiR() field
             as XiR()[0], XiR()[1], XiR()[2]

    phi itself is a particle-carried scalar that is mixed pairwise by the
    first-conditioning model and relaxed towards the Eulerian progress
    variable c by MixingPopeCloud::updatePhi() over a timescale tauC, so it
    stays anchored to the resolved solution instead of free-running.

    Pairing uses a k-d tree search local to this class (findPairsSC /
    kdTreeSearchSC) in which the three PHYSICAL coordinates, normalised by
    r_i, compete with the four reference-space axes, normalised by the Xim_i
    entries.  This mirrors premixedMixParticleModel rather than the inherited
    mixParticleModel::KkdTreeLikeSearch, whose physical-space branches are
    commented out.

    That distinction matters.  In MMCcurl the reference variables ARE the
    shadow positions, so they supply locality of mixing themselves.  Here the
    reference space mixes phi_modified against the shadow positions, so if a
    tight phiMod_m makes phi_modified dominant there is nothing left to keep
    pairs local in physical space, and the model degenerates into mixing
    composition between particles at opposite ends of the domain.  Keeping the
    physical axes in contention means phiMod_m can be tightened to give the
    progress variable more weight WITHOUT sacrificing locality.

    Because the normalisers divide, a SMALLER value makes an axis MORE likely
    to be chosen as the split dimension.

    The aISO timescale is computed identically to MMCcurl (from CE only; CL is
    accepted for dictionary compatibility but unused).  Only Y, T, and hA are
    mixed between the paired particles; phi, XiR, XiC, and secondCondFlag are
    NOT mixed.

    Parallel pairing is not implemented; local pairing only.

    Diagnostics are written to <case>/postProcessing/secondCondMMCcurl.log:
    flagged-particle counts, pair/triple counts, the split-axis histogram over
    all seven candidate axes, and the realised phiModified range.  Note that
    phiModified is NOT bounded by 1 even though phi is, which matters when the
    thermophysical coupling conditions on it (see KernelEstimation's fHigh).

    Usage in cloudProperties:
    \verbatim
    mixingModel    secondCondMMCcurl;

    secondCondMMCcurlCoeffs
    {
        // Physical-space normalisation. Live: it sets how strongly the k-d
        // tree insists on locality of mixing against the reference axes.
        r_i     1e-3;

        // 4-entry normalisation for the (phiModified, xi_x, xi_y, xi_z)
        // reference space.  Tight phiMod_m gives the progress variable more
        // weight in the splitting.
        Xim_i
        {
            phiMod_m    0.01;
            sPx_m       0.05;
            sPy_m       0.05;
            sPz_m       0.05;
        }

        pairingMethod   local;

        // aISO timescale coefficients
        CE              0.1;
        meanTimeScale   true;
    }
    \endverbatim

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
 
 
    // Private member functions
 
        //- Mix only Y, T, and hA between paired particles.
        //  Deliberately does NOT mix phi, XiR, XiC, or secondCondFlag.
        static void mixSpeciesOnly
        (
            particleType& p,
            particleType& q,
            scalar mixExtent
        );

        //- Pair the flagged particles using a k-d tree search in which the
        //  three PHYSICAL coordinates compete against the four reference-space
        //  axes for the split.
        //
        //  The shared mixParticleModel::findPairs()/KkdTreeLikeSearch() has its
        //  physical-space branches commented out, so only XiR axes are split
        //  candidates. That is correct for MMCcurl, where the reference
        //  variables ARE the shadow positions and therefore supply locality
        //  themselves. It is not correct here: the second-conditioning space
        //  mixes phi_deg against the shadow positions, so any normalisation
        //  that makes phi_deg dominant removes the only locality mechanism and
        //  the model starts pairing particles that are arbitrarily far apart
        //  in physical space.
        //
        //  These two routines mirror premixedMixParticleModel, which keeps
        //  disX/disY/disZ (normalised by r_i) as ncond 0,1,2 alongside the
        //  reference variables at ncond = 3+i. With them, tightening phiMod_m
        //  raises the weight given to phi_deg WITHOUT ever sacrificing
        //  physical locality, and r_i in secondCondMMCcurlCoeffs becomes live.
        //
        //  The shared base implementation is left untouched.
        void findPairsSC
        (
            const DynamicList<eulerianFieldData>& eulerianFieldList,
            DynamicList<List<label>>& pairs
        ) const;

        void kdTreeSearchSC
        (
            const DynamicList<eulerianFieldData>& particleList,
            label l,
            label u,
            std::vector<label>& pInd,
            std::vector<label>& L,
            std::vector<label>& U
        ) const;
 
 
protected:
 
    // Override of base-class particle-list builder
 
        //- Build the local particle list filtered to secondCondFlag == 1,
        //  populate eulerianFieldDataList_ with a 4D XiR
        //  (phiModified, xi_x, xi_y, xi_z), interpolate Eulerian fields,
        //  and call findPairs() for local pairing.
        virtual void buildParticleList() override;
 
 
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
 
        //- Boundary initialisation (not applicable — returns empty field)
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

