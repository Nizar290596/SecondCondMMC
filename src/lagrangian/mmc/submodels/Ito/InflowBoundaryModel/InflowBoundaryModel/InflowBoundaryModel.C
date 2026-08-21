/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "InflowBoundaryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CloudType>
Foam::InflowBoundaryModel<CloudType>::InflowBoundaryModel(CloudType& owner)
    :
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::InflowBoundaryModel<CloudType>::InflowBoundaryModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
    :
    SubModelBase<CloudType>(owner, dict, typeName,type)
{}


template<class CloudType>
Foam::InflowBoundaryModel<CloudType>::InflowBoundaryModel
(
    const InflowBoundaryModel& ibc
)
    :
    SubModelBase<CloudType>(ibc)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InflowBoundaryModel<CloudType>::~InflowBoundaryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InflowBoundaryModel<CloudType>::inflow()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InflowBoundaryModelNew.C"

// ************************************************************************* //

