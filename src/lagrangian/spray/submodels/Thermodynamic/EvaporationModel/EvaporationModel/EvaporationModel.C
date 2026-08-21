/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "EvaporationModel.H"
#include "mathematicalConstants.H"
#include "meshTools.H"
#include "volFields.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EvaporationModel<CloudType>::EvaporationModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner)

{}
    


template<class CloudType>
Foam::EvaporationModel<CloudType>::EvaporationModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type)
{
}


template<class CloudType>
Foam::EvaporationModel<CloudType>::EvaporationModel
(
    const EvaporationModel<CloudType>& evm
)
:
    CloudSubModelBase<CloudType>(evm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EvaporationModel<CloudType>::~EvaporationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
Foam::scalarField Foam::EvaporationModel<CloudType>::YFuelVap
(
     const scalarField YFuel,
     const scalar Ts,
     const scalar Pa,
     const scalar Lv,
     const scalar molWtSstate
) const
{
    notImplemented
    (
        "Foam::scalar Foam::EvaporationModel<CloudType>::fSurface"
        "("
            "const scalarField,"
            "const scalar,"
            "const scalar,"
            "const scalar,"
            "const scalar"
        ") const"
    );

    scalarField YDummy;
    YDummy.setSize(YFuel.size(),0.0);
    return YDummy;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EvaporationModelNew.C"


// ************************************************************************* //

