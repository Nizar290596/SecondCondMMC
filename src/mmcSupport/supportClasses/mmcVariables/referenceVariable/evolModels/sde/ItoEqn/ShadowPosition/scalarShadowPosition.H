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

Typedef
    Foam::scalarShadowPosition
Description
    A scalar version of the templated ShadowPosition
SourceFiles
    scalarShadowPosition.C
\*---------------------------------------------------------------------------*/

#ifndef scalarShadowPosition_H
#define scalarShadowPosition_H

#include "scalar.H"
#include "ItoEqn.H"
#include "scalarItoEqn.H"
#include "interpolation.H"
#include "HashTable.H"
//#include "ShadowPosition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class T>
class ShadowPosition;

template<>
class ShadowPosition<scalar>;

typedef ShadowPosition<scalar> scalarShadowPosition;

template<class T>
scalar LESomegaSP
(
    const vector& position,
    //const vector& coordinates,
    label cell,
    label face,
    ShadowPosition<T>& sp
);

template<class T>
scalar RASomegaSP
(
    const vector& position,
    //const vector& coordinates,
    label cell,
    label face,
    ShadowPosition<T>& sp
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\
                 Class ShadowPosition Declaration
\*---------------------------------------------------------------------------*/

template<>
class ShadowPosition<scalar>
:
    public scalarItoEqn
{
    public:

        typedef scalar (*timeScaleModelPtr)
            (
                const vector& position,
		//const vector& coordinates,
                label cell,
                label face,
                const scalarShadowPosition& sp
            );

//        typedef HashTable<timeScaleModelPtr> timeScaleModelTable;

//        timeScaleModelTable timeScaleModelTable_;

        static HashTable<timeScaleModelPtr> timeScaleModelTable_;

    protected:

        // Protected data

            //- Drift Coefficient constant 
            scalar a_;

            //- Diffusion Coefficient constant
            scalar b_;

            //- Variable lambda calculation
            Switch variableLambda_;

            //- shadow component (x or y or z)
            word componentName_;

            //- shadow component (x or y or z)
            label component_;

            //- Turbulence model type (LES or RAS or laminar)
            word modelType_;

            //- Interpolation schemes
            dictionary intpSchemes_;

//            autoPtr<interpolation<scalar>> UInterp_;
            autoPtr<interpolation<vector>> UInterp_;

            autoPtr<interpolation<scalar>> DtInterp_;

            autoPtr<interpolation<scalar>> DInterp_;

            autoPtr<interpolation<scalar>> DeltaEInterp_;
            
            autoPtr<interpolation<scalar>> vbInterp_;  


    public:

        //- Runtime type information
        TypeName("scalarShadowPosition");

        // Constructors

            //- Construct from dictionary and mesh
            ShadowPosition
            (
                const dictionary& entry,
                const fvMesh& mesh
            );

            //- Construct copy
            ShadowPosition
            (
                const ShadowPosition& sp
            );

        //- Destructor
        virtual ~ShadowPosition()
        {};

        // Member functions

            static void constructTimeScaleModelTable(); //setModels();

            // Access

                //- Return Drift Coefficient
                inline const scalar& a() const;

                //- Return Diffusion Coefficient
                inline const scalar& b() const;

            // Edit

                //- Return Drift Coefficient
                inline scalar a();

                //- Return Diffusion Coefficient
                inline scalar b();

            //- computation of coefficients

                friend scalar LESomegaSP
                (
                    const vector& position,
		    //const vector& coordinates,
                    label cell,
                    label face,
                    const scalarShadowPosition& sp
                );

                friend scalar RASomegaSP
                (
                    const vector& position,
		    //const vector& coordinates,
                    label cell,
                    label face,
                    const scalarShadowPosition& sp
                );

                scalar omegaSPos
                (
                    const vector& position, 
		    //const vector& coordinates,
                    label cell, 
                    label face
                ) const;

                virtual void computeDriftCoeff
                (
                    const particle& p,
                    const scalar& x
                ) const;

                virtual void computeDiffusionCoeff(const particle& p) const;

            //- Intial conditions
                virtual scalar t0
                (
                    const particle& p,
                    scalar dt, 
                    const scalar& x
                ) const;

                virtual autoPtr<ItoEqn<scalar>> clone() const
                {
                    return autoPtr<ItoEqn<scalar>>(new ShadowPosition<scalar>(*this));
                }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
