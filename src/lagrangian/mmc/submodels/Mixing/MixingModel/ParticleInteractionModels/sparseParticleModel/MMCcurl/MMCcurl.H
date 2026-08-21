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
    Foam::MMCcurl

Description
    MMC mixing based on Curl's model with a K-D tree like mixing pair structure
    
    Detail can be found in:
    [1] Son Vo, O.T. Stein, A. Kronenburg, M.J. Cleary,
        Assessment of mixing time scales for a sparse particle method,
        Combustion and Flame,
        Volume 179,
        2017,
        Pages 280-299,
        ISSN 0010-2180,
        https://doi.org/10.1016/j.combustflame.2017.02.017. 
    

\*---------------------------------------------------------------------------*/

#ifndef MMCcurl_H
#define MMCcurl_H

#include "mixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MMCcurl
:
    public mixParticleModel<CloudType>
{
    // Typedef for particle type
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;
    
    // Private data
    
        //- Time scale model coefficients
        const scalar CL_;
        const scalar CE_;
        const scalar beta_;

        //- aISO control
        const Switch aISO_;

        //- true: harmonic mean of two particle mixing time scale
        //- false: min of two particle mixing time scale
        const Switch meanTimeScale_;
    
    public:

        //- Runtime type information
        TypeName("MMCcurl");

    // Constructors

        //- Construct from dictionary
        MMCcurl(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        MMCcurl(const MMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new MMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~MMCcurl() = default;
        

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
        
        virtual const scalarField XiR0(label patch, label patchFace, particle& p);
        virtual const scalarField XiR0(label celli, particle& p);   

    // Info Statement

        //- Print an info statement to the screen
        void printInfo();   

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
