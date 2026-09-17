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
    Sparse species mixing conditioned strongly on phiModified and weakly on
    physical position. Optional shadow coordinates are resolved by name.
    includeShadowPositions selects (phiModified,x) or (phiModified,xi,x).
    Only particles with secondCondFlag == 1 participate. Species, absolute
    enthalpy and coupling scalars are mixed conservatively; temperature is
    reconstructed. The dense reference and selection flag are not mixed here.

    See docs/secondConditioning.md for normalization, rejection, random mixing,
    prescribed/aISO timescales, parallel pairing, and validation requirements.

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
        const scalar CE_;
        const Switch includeShadowPositions_;
        const label diagnosticInterval_;
        void configureReferenceSpace();
 
        //- true  → harmonic mean of the two particle mixing time scales
        //- false → minimum of the two particle mixing time scales
        const Switch meanTimeScale_;
 
 
    // Private member functions
 
        //- Mix Y, hA and XiC; reconstruct T from thermodynamics.
        //  Does not mix phi, XiR, burnedAge or secondCondFlag.
        void mixSpeciesOnly
        (
            particleType& p,
            particleType& q,
            scalar mixExtent
        );
 
 
protected:
 
    // Override of base-class particle-list builder
 
        //- Build the local particle list filtered to secondCondFlag == 1,
        //  populate the selected reference coordinates and transport fields,
        //  then use the selected local/global/neighbour pairing method.
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

