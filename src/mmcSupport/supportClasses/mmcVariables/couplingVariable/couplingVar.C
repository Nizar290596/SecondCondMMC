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

#include "couplingVar.H"

namespace Foam
{
    defineTypeNameAndDebug(couplingVar, 0);
    addToRunTimeSelectionTable(basicMMCVar,couplingVar,dictionary);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//inline Foam::word Foam::couplingVar::cName()
//{
//    return cName_;
//}               

//inline Foam::Switch Foam::couplingVar::passive()
//{
//    return passive_;
//}            


//inline Foam::label Foam::couplingVar::mixStep()
//{
//    return mixStep_;
//}  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::couplingVar::couplingVar
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    basicMMCVar(entry, mesh),

    cName_(entry.lookup("couplingName")),
    
    passive_(entry.lookup("passive")),
    
    mixStep_(readLabel(entry.lookup("mixStep")))
{}


Foam::couplingVar::couplingVar(const couplingVar& cm)
:
    basicMMCVar(cm),

    cName_(cm.cName_),

    passive_(cm.passive_),

    mixStep_(cm.mixStep_)
{}
// ************************************************************************* //
