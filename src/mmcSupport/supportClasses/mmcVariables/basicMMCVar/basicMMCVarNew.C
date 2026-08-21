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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicMMCVar>
Foam::basicMMCVar::New
(
    const dictionary& entry,
    const fvMesh& mesh
)
{
    const word variableType(entry.lookup("MMCType"));

//    Info << "Creating a mmc Variable named: "<< word(entry.lookup("name")) 
//         << " of type: " << variableType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(variableType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicMMCVar::New"
            "("
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Unknown mmc variable type "
            << variableType << nl << nl
            << "Valid mmc variable types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicMMCVar>(cstrIter()(entry, mesh));
}


// ************************************************************************* //
