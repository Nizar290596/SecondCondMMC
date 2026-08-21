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

#include "referenceVar.H"

namespace Foam
{
    defineTypeNameAndDebug(referenceVar, 0);
    addToRunTimeSelectionTable(basicMMCVar,referenceVar,dictionary);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//inline Foam::word Foam::referenceVar::rName()
//{
//    return rName_;
//}               

//inline Foam::word Foam::referenceVar::refType()
//{
//    return refType_;
//}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::referenceVar::referenceVar
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    basicMMCVar(entry, mesh),
    
    rName_(this->name()),// reference name same as LES field

    refType_(entry.lookup("referenceType")),

    evolMethod_
    (
	evolModel::New
	(
	    refType_,
	    rName_,
	    mesh
	)
    )
{}

Foam::referenceVar::referenceVar(const referenceVar& cm)
:
    basicMMCVar(cm),

    rName_(cm.rName_),

    refType_(cm.refType_),

    evolMethod_(cm.evolMethod_().clone())
{}
// ************************************************************************* //
