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

#include "SynthesisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SynthesisModel<CloudType>::SynthesisModel(CloudType& owner)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::SynthesisModel<CloudType>::SynthesisModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type)
{}


template<class CloudType>
Foam::SynthesisModel<CloudType>::SynthesisModel
(
    const SynthesisModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SynthesisModel<CloudType>::~SynthesisModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SynthesisModel<CloudType>::synthesize()
{}


template<class CloudType> 
const Foam::wordList Foam::SynthesisModel<CloudType>::psdPropertyNames() const
{
    wordList psdPropertyNames;
    psdPropertyNames.setSize(0);
    
    return psdPropertyNames;
}

template<class CloudType> 
const Foam::wordList Foam::SynthesisModel<CloudType>::physicalAerosolPropertyNames() const
{
    wordList physicalAerosolPropertyNames;
    physicalAerosolPropertyNames.setSize(0);
    
    return physicalAerosolPropertyNames;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SynthesisModelNew.C"

// ************************************************************************* //

