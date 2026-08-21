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

\*---------------------------------------------------------------------------*/

#include "basicMMCVar.H"

namespace Foam
{
    defineTypeNameAndDebug(basicMMCVar, 0);
    defineRunTimeSelectionTable(basicMMCVar, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMMCVar::basicMMCVar
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    name_(entry.lookup("name")),

    field_(nullptr)
{
    IOobject header
    (
        name_,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );
    
    if (header.typeHeaderOk<volScalarField>(true))
    {
        field_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info << "File "<<name_<<" could not be found in "
             << mesh.time().timeName() << nl
             << "Using Xidefault in " << mesh.time().timeName() << endl;

        volScalarField Xidefault
        (
            IOobject
            (
                "Xidefault",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        field_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Xidefault
            )
        );
    }
}

Foam::basicMMCVar::basicMMCVar(const basicMMCVar& bmv)
:
    name_(bmv.name_),
    field_(new volScalarField(bmv.field_()))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicMMCVarNew.C"

// ************************************************************************* //
