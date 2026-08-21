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

#include "CloudMixingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudMixingModel<CloudType>::CloudMixingModel(CloudType& owner)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::CloudMixingModel<CloudType>::CloudMixingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type)
{}


template<class CloudType>
Foam::CloudMixingModel<CloudType>::CloudMixingModel
(
    const CloudMixingModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudMixingModel<CloudType>::Smix
(
) 
{
    NotImplemented;
}


template<class CloudType>
inline const Foam::scalarField Foam::CloudMixingModel<CloudType>::XiR0(label patch, label patchFace, particle& p) 
{
    NotImplemented;
    return scalarField::null();
}


template<class CloudType>
inline const Foam::scalarField Foam::CloudMixingModel<CloudType>::XiR0(label celli, particle& p) 
{
    NotImplemented;
    return scalarField::null();
}


template<class CloudType>
inline const Foam::wordList& Foam::CloudMixingModel<CloudType>::XiRNames() const
{
    NotImplemented;
    return wordList::null();//!!!
}


template<class CloudType>
inline Foam::label Foam::CloudMixingModel<CloudType>::numXiR() const
{
    NotImplemented;
    return 0;
};


template<class CloudType>
inline const Foam::mmcVarSet& Foam::CloudMixingModel<CloudType>::XiR() const
{
    NotImplemented;
    return mmcVarSet::null();
}

// ************************************************************************* //

