/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "ThermoPhysicalCouplingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::ThermoPhysicalCouplingModel<CloudType> >
Foam::ThermoPhysicalCouplingModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
{
    const word modelType(dict.lookup("thermoPhysicalCouplingModel"));

    Info << nl << "Selecting density coupling model " << modelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    //    Info <<"I'm here"<<endl;
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ThermoPhysicalCouplingModel<CloudType>::New"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Unknown density coupling model type "
            << modelType << nl << nl
            << "Valid density coupling model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<ThermoPhysicalCouplingModel<CloudType> >(cstrIter()(dict, owner, Xi));
}


// ************************************************************************* //
