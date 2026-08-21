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
    Foam::ItoEqn

Description
    Generic version of a Ito SDE

SourceFiles
    ItoEqn.C

\*---------------------------------------------------------------------------*/

#ifndef ItoEqn_H
#define ItoEqn_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"

#include "StochasticLib.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class ItoEqn Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class ItoEqn
{
    typedef T ItoType;

    protected:

        // Protected data

            //- Drift Coefficient
            mutable T A_;

            //- Diffusion Coefficient
            mutable T D_;

            //- Random numbers
            mutable StochasticLib1 rndGen_;

    public:

        //- Runtime type information
        TypeName("ItoEqn");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ItoEqn,
            dictionary,
            (
                const dictionary& entry,
                const fvMesh& mesh
            ),
            (entry, mesh)
        );

        // Constructors

            //- Null Construct
            ItoEqn();


        //- Destructor
        virtual ~ItoEqn(){};

        //- Selector
        static autoPtr<ItoEqn> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        // Member functions

            // Access

                //- Return Drift Coefficient
                inline const T& A() const;

                //- Return Diffusion Coefficient
                inline const T& D() const;

            // Edit

                //- Return Drift Coefficient
                inline T& A();

                //- Return Diffusion Coefficient
                inline T& D();

            //- virtual functions
                virtual void computeDriftCoeff
                (
                    const particle& p, 
                    const T& x
                ) const;

                virtual void computeDiffusionCoeff(const particle& p) const;

            //- Integration
                virtual T integrate
                (
                    const particle& p,
                    scalar dt, 
                    const T& x
                ) const;

            //- Intial conditions
                virtual T t0
                (
                    const particle& p,
                    scalar dt, 
                    const T& x
                ) const = 0;

//                virtual T t0
//                (
//                    const particle& p,
//                    scalar dt, 
//                    const T& x
//                ) const = 0;

            //- Construct and return a clone to Ito Eqn
                virtual autoPtr<ItoEqn<T>> clone() const = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeItoEqn(ItoVariableType)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(ItoEqn<ItoVariableType>,0);            \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        ItoEqn<ItoVariableType>,                                               \
        dictionary                                                             \
    );

#define makeItoEqnType(SS,ItoVariableType)                                     \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<ItoVariableType>,0);                \
                                                                               \
    ItoEqn<ItoVariableType>::                                                  \
        adddictionaryConstructorToTable<SS<ItoVariableType> >                  \
            add##SS##ItoVariableType##ConstructorToTable_;

#define makeSpecializedItoEqnType(SS,ItoVariableType)                          \
                                                                               \
    defineTypeNameWithName                                                     \
    (                                                                          \
        SS<ItoVariableType>,                                                   \
        SS<ItoVariableType>::typeName_()                                       \
    );                                                                         \
                                                                               \
    defineDebugSwitchWithName                                                  \
    (                                                                          \
        SS<ItoVariableType>,                                                   \
        SS<ItoVariableType>::typeName_(),                                      \
        0                                                                      \
    );                                                                         \
                                                                               \
    registerTemplateDebugSwitchWithName                                        \
    (                                                                          \
        SS<ItoVariableType>,                                                   \
        SS<ItoVariableType>::typeName_()                                       \
    );                                                                         \
                                                                               \
    ItoEqn<ItoVariableType>::                                                  \
        adddictionaryConstructorToTable<SS<ItoVariableType> >                  \
            add##SS##ItoVariableType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ItoEqnI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ItoEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
