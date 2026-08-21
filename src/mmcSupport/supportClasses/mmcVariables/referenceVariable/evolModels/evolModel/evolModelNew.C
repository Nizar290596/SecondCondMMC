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

#include "evolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::evolModel> Foam::evolModel::New
(
    const word& refType,
    const word& refName,
    const fvMesh& mesh
)
{
//    typename wordConstructorTable::iterator cstrIter =
//        wordConstructorTablePtr_->find(refType);

    auto cstrIter =
        wordConstructorTablePtr_->find(refType);


    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "evolModel::New"
            "("
                "const word&"
                "const word&"
                "const fvMesh&"
            ")"
        )   << "Unknown evolution model for reference variables (refType)"
            << refType << nl << nl
            << "Valid refType are:" << nl
            << wordConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<evolModel>(cstrIter()(refName, mesh));
}

// ************************************************************************* //
