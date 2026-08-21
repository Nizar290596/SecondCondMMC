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
      - phi_modified = phiModified (modified progress variable, already
                       computed from the OU process at each timestep)
      - xi = shadow positions stored in the particle's existing XiR() field
             as XiR()[0], XiR()[1], XiR()[2]
 
    The k-d tree splitting uses the normalisation factors from the inherited
    Xii_ list (read from the Xim_i subdictionary — must have exactly 4
    entries corresponding to phiMod_m, sPx_m, sPy_m, sPz_m).  A tight
    normalisation for phi_modified (e.g. 0.01) causes the tree to split
    preferentially on the progress-variable axis; looser normalisations for
    the shadow-position coordinates (e.g. 0.05) provide spatial context with
    lower priority.  Physical coordinates are NOT used as split dimensions
    (the ncond = 3 + i offset in KkdTreeLikeSearch ensures this).
 
    The aISO timescale is computed identically to MMCcurl.  Only Y, T, and
    hA are mixed between the paired particles; phi, XiR, XiC, and
    secondCondFlag are NOT mixed.
 
    Parallel pairing is not implemented; local pairing only.
 
    Usage in cloudProperties:
    \verbatim
    mixingModel    secondCondMMCcurl;
 
    secondCondMMCcurlCoeffs
    {
        // Base class: position-space normalisation (unused for k-d splitting
        // but required by the constructor)
        r_i     1e-3;
 
        // 4-entry normalisation for the (phiModified, xi_x, xi_y, xi_z)
        // reference space.  Tight phiMod_m dominates the k-d tree splitting.
        Xim_i
        {
            phiMod_m    0.01;
            sPx_m       0.05;
            sPy_m       0.05;
            sPz_m       0.05;
        }
 
        pairingMethod   local;
 
        // aISO timescale coefficients (same meaning as MMCcurl)
        CL              0.5;
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

