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
    Foam::sde

Description
    This class is used for the evolution of the referece variable using SDE
    in the Ito sense.

SourceFiles
    sde.C

\*---------------------------------------------------------------------------*/

#ifndef sde_H
#define sde_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"

#include "evolModel.H"
#include "ItoEqn.H"
#include "scalarItoEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sde Declaration
\*---------------------------------------------------------------------------*/

class sde
:
    public evolModel
{
    protected:

        // Protected data

        autoPtr<ItoEqn<scalar> > SDEModel_;

    public:

        //- Runtime type information
        TypeName("sde");

        // Constructors

            //- Construct from dictionary and mesh
            sde
            (
                const word& refName,
                const fvMesh& mesh
            );

            //- Copy Constructor
            sde(const sde& sd);


        //- Destructor
        virtual ~sde()
        {};

        // Member functions

            //- virtual functions
                virtual scalar compute
                (
                    const particle& p,
                    scalar dt, 
                    const scalar& x
                ) const;

                virtual scalar initialise
                (
                    label patch,
                    label patchFace,
                    particle& p
                ) const;

                virtual scalar initialise
                (
                    label celli,
                    particle& p
                ) const;

            //- Construct and return a clone to evolModel
                virtual autoPtr<evolModel> clone() const
                {
                    return autoPtr<evolModel>(new sde(*this));
                }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
