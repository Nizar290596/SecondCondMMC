/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ShadowPosition

Description
    Generic version of a ShadowPosition SDE. 

    Reference: A model for turbulent mixing based on shadow-position 
               conditioning.

SourceFiles
    ShadowPosition.C

\*---------------------------------------------------------------------------*/

#ifndef ShadowPosition_H
#define ShadowPosition_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "ItoEqn.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class T>
class ShadowPosition;

template<>
class ShadowPosition<scalar>;

/*---------------------------------------------------------------------------*\
                 Class ShadowPosition Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class ShadowPosition
:
    public ItoEqn<T>
{
    protected:

        // Protected data

            //- Drift Coefficient constant 
            scalar a_;

            //- Diffusion Coefficient constant
            scalar b_;
            
            //- Variable lambda calculation
            Switch variableLambda_;
            
            //- Turbulence model type (LES or RANS)
            word modelType_;

            autoPtr<interpolation<T>> UInterp_;

            autoPtr<interpolation<scalar>> DEffInterp_;

            autoPtr<interpolation<scalar>> DeltaEInterp_;

            autoPtr<interpolation<scalar>> vbInterp_;            

    public:

        //- Runtime type information
        TypeName("ShadowPosition");

        // Constructors

            //- Construct from dictionary and mesh
            ShadowPosition
            (
                const dictionary& entry,
                const fvMesh& mesh
            );

        //- Destructor
        virtual ~ShadowPosition(){};

        // Member functions

            // Access

                //- Return Drift Coefficient
                inline scalar a() const;

                //- Return Diffusion Coefficient
                inline scalar b() const;

            // Edit

                //- Return Drift Coefficient
                inline scalar& a();

                //- Return Diffusion Coefficient
                inline scalar& b();

            //- virtual functions
//                template<class ParticleType>
                void computeDriftCoeff
                (
                    const T& x,
                    const particle& p
                ) const;

//                template<class ParticleType>
                void computeDiffusionCoeff(const particle& p) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ShadowPositionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ShadowPosition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
