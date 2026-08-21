/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-217 OpenFOAM Foundation
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

#include "crVar.H"

namespace Foam
{
    defineTypeNameAndDebug(crVar, 0);
    addToRunTimeSelectionTable(basicMMCVar,crVar,dictionary);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::crVar::crVar
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    basicMMCVar(entry, mesh),

  //  couplingVar(entry, mesh),
    
    referenceVar(entry, mesh),

    couplingVar(entry, mesh)
{}

Foam::crVar::crVar
(
    const crVar& cm
)
:
    basicMMCVar(cm),

    referenceVar(cm),

    couplingVar(cm)
{}
// ************************************************************************* //
