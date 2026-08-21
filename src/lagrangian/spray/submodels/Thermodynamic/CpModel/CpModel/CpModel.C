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

#include "CpModel.H"
#include "mathematicalConstants.H"
#include "meshTools.H"
#include "volFields.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CpModel<CloudType>::CpModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::CpModel<CloudType>::CpModel
(
    const dictionary& dict,
    CloudType& owner, 
    const word& modelType
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, modelType)
{
    Info<< "\n Constructing Cp model base!" << endl;
}


template<class CloudType>
Foam::CpModel<CloudType>::CpModel
(
    const CpModel<CloudType>& cpm
)
:
    CloudSubModelBase<CloudType>(cpm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CpModel<CloudType>::~CpModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CpModelNew.C"

// ************************************************************************* //
