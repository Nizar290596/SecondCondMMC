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

#include "particleMatchingAlgorithm.H"
#include "kdTreeLikeMatching.H"
#include "greedyBiPartiteMatching.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class particleType>
Foam::autoPtr<Foam::particleMatchingAlgorithm<particleType>>
Foam::particleMatchingAlgorithm<particleType>::New
//Foam::particleMatchingAlgorithm<eulerianFieldData>::New
(
    const scalar& ri, 
    const List<scalar>& Xii,
    const dictionary& dict
)
{
    const word modelType(dict.getOrDefault<word>("particleMatchingAlgorithm","kdTreeLike"));

    Info << "Selecting particle matching algorihtm " << modelType << endl;

    autoPtr<particleMatchingAlgorithm<particleType>> model;
    //autoPtr<particleMatchingAlgorithm<eulerianFieldData>> model

    if (modelType == "kdTreeLike")
    {
        model.reset
        (
            new kdTreeLikeMatching<particleType>(ri,Xii,dict)
            //new kdTreeLikeMatching<eulerianFieldData>(ri,Xii,dict)
        );

    }
    else if (modelType == "greedyBiPartiteMatching")
    {
        model.reset
        (
            new greedyBiPartiteMatching<particleType>(ri,Xii,dict)
            //new greedyBiPartiteMatching<eulerianFieldData>(ri,Xii,dict)
        );
    }
    else
    {
        FatalErrorIn
        (
            "particleMatchingAlgorithm::New"
            "()"
        )   << "Unknown mixing model type "
            << modelType << nl << nl
            << "Valid mixing model types are:" << nl
            << "\tkdTreeLike" << nl
            << "\tgreedyBiPartiteMatching" << nl
            << exit(FatalError);
    }

    return model;
}


// ************************************************************************* //
