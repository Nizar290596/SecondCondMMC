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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::autoPtr<Foam::ItoEqn<T>>
Foam::ItoEqn<T>::New
(
    const dictionary& entry,
    const fvMesh& mesh
)
{
    const word SDEType(entry.lookup("SDEModel"));

//    typename dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(SDEType);

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(SDEType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ItoEqn::New"
            "("
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Unknown ItoEqn type "
            << SDEType << nl << nl
            << "Valid ItoEqn types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<ItoEqn<T>>(cstrIter()(entry, mesh));
}

// ************************************************************************* //
