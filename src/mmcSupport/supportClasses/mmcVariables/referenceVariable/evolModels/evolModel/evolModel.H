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
    Foam::evolModel

Description
    This class is the base/interface to handle reference variables evolution
    models.

SourceFiles
    evolModel.C

\*---------------------------------------------------------------------------*/

#ifndef evolModel_H
#define evolModel_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"

#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class evolModel Declaration
\*---------------------------------------------------------------------------*/

class evolModel
{
    protected:

        // Protected data
        IOdictionary referenceVariableProperties_;

    public:

        //- Runtime type information
        TypeName("evolModel");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            evolModel,
            word,
            (
                const word& refName,
                const fvMesh& mesh
            ),
            (refName, mesh)
        );

        // Constructors

            //- Construct from dictionary and mesh
            evolModel(const fvMesh& mesh);

            //- Copy Constructor
            evolModel(const evolModel& em);


        //- Destructor
        virtual ~evolModel()
        {};

        //- Selector
        static autoPtr<evolModel> New
        (
            const word& refType,
            const word& refName,
            const fvMesh& mesh
        );

        // Member functions

            //- Edit
                IOdictionary& rVarProperties()
                {
                    return referenceVariableProperties_;
                }

            //- Access
                const IOdictionary& rVarProperties() const
                {
                    return referenceVariableProperties_;
                }

            //- virtual functions

                virtual scalar compute
                (
                    const particle& p,
                    scalar dt,
                    const scalar& x
                ) const = 0;

                virtual scalar initialise
                (
                    label patch,
                    label patchFace,
                    particle& p
                ) const = 0;

                virtual scalar initialise
                (
                    label celli,
                    particle& p
                ) const = 0;

            //- Null evolution model
                static evolModel& null()
                {
                    return *reinterpret_cast<evolModel* >(0);
                }

            //- Construct and return a clone
            virtual autoPtr<evolModel> clone() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
