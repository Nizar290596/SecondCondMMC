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

    Pairing uses the inherited mixParticleModel::findPairs() /
    KkdTreeLikeSearch().  Splitting is therefore over the FOUR reference axes
    only — phiModified and the three shadow positions — and NEVER over the
    physical coordinates: the physical-space branches in KkdTreeLikeSearch are
    deliberately disabled, and the ncond = 3 + i offset points lessArg past
    x, y and z.  This matches MMCcurl, where the shadow positions are likewise
    the sole carriers of locality.

    Consequently r_i is NOT a split axis here.  It is still a mandatory
    dictionary entry because the base-class constructor reads it (and passes it
    to mixingSubVolumes), but with pairingMethod local it takes no part in any
    pairing decision.

    The balance between conditioning on reaction progress and remaining local
    in shadow space is therefore set ENTIRELY by the ratio of the Xim_i
    entries: phiMod_m against sPx_m / sPy_m / sPz_m.  Because the normalisers
    divide, a SMALLER value makes an axis MORE likely to be chosen.

    Note that the k-d tree compares ROOT extents, not local gradients: at the
    top node the shadow axes span the whole cloud while phiModified spans its
    own (order-unity) range.  With three shadow axes against one phiModified
    axis, equal root normalised extents give phiModified about a quarter of the
    splits.  Setting sP*_m here to the values validated for MMCcurl will make
    the shadow axes utterly dominant and phiModified will essentially never be
    selected — these Xim_i are read from this model's own dictionary and are
    meant to be looser than the first-conditioning ones.

    The aISO timescale is computed identically to MMCcurl (from CE only; CL is
    accepted for dictionary compatibility but unused).  Only Y, T, and hA are
    mixed between the paired particles; phi, XiR, XiC, and secondCondFlag are
    NOT mixed.

    Parallel pairing is not implemented; local pairing only.

    Diagnostics are written to <case>/postProcessing/secondCondMMCcurl.log:
    flagged-particle counts, pair/triple counts, the split-axis histogram over
    the four reference axes, and the realised phiModified range.  Note that
    phiModified is NOT bounded by 1 even though phi is, which matters when the
    thermophysical coupling conditions on it (see KernelEstimation's fHigh).

    Usage in cloudProperties:
    \verbatim
    mixingModel    secondCondMMCcurl;

    secondCondMMCcurlCoeffs
    {
        // Required by the base-class constructor, but NOT a split axis:
        // physical coordinates never take part in the pairing.
        r_i     1e-3;

        // 4-entry normalisation for the (phiModified, xi_x, xi_y, xi_z)
        // reference space -- the only thing that sets the balance between
        // conditioning on progress and staying local in shadow space.
        // These are this model's OWN values: they are read from
        // secondCondMMCcurlCoeffs, not from MMCcurlCoeffs, and are meant to
        // be looser than the first-conditioning ones.
        Xim_i
        {
            phiMod_m    0.1;
            sPx_m       1.3e-3;
            sPy_m       1.3e-3;
            sPz_m       1.3e-3;
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

