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
    Foam::premixedMMCcurl

Description
    premixed MMC mixing based on Curl's model with a K-D tree like mixing pair structure
    
    Detail can be found in:
    [1] Iaroslavtceva, N., Kronenburg, A. & Stein, O.T,
        Multiple Mapping Conditioning Mixing Time Scales for Turbulent Premixed Flames,
        Flow Turbulence Combust,
        Volume 110,
        2023,
        Pages 395–415,
        https://doi.org/10.1007/s10494-022-00375-1. 
    

\*---------------------------------------------------------------------------*/

#ifndef premixedMMCcurl_H
#define premixedMMCcurl_H

#include "premixedMixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class premixedMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class premixedMMCcurl
:
    public premixedMixParticleModel<CloudType>
{
    // Typedef for particle type
    typedef typename CloudType::particleType particleType;
    typedef typename premixedMixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;
    
    // Private data
    
        //- Mixing time scale model
        const word mixTimeScale_;
    
        //- Time scale model coefficients
        const scalar CL_;
        const scalar CE_;
        const scalar Ka_;
        const scalar alpha_;

        //- Limit mixing of particles that are far in composition space 
        //- true: particles with dXiR > 2*Ximi are NOT mixed
        //- false: no limiting of mixing based on dXiR
        const Switch limitkdTree_;

        //- true: harmonic mean of two particles mixing time scales
        //- false: min of two particles mixing time scales
        const Switch meanTimeScale_;
    
    public:

        //- Runtime type information
        TypeName("premixedMMCcurl");

    // Constructors

        //- Construct from dictionary
        premixedMMCcurl(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        premixedMMCcurl(const premixedMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new premixedMMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~premixedMMCcurl() = default;
        

    // Member Functions
    
        //- Mixing particle pairs
        //  See also [1] for description of the time scale
        void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        );
        
        virtual const scalarField XiR0(label patch, label patchFace);
        virtual const scalarField XiR0(label celli);   

    // Info Statement

        //- Print an info statement to the screen
        void printInfo();   

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "premixedMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
