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

#include "SootModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SootModel<CloudType>::SootModel(CloudType& owner)
:
    SubModelBase<CloudType>(owner)
{}

template<class CloudType>
Foam::SootModel<CloudType>::SootModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type)
{}

template<class CloudType>
Foam::SootModel<CloudType>::SootModel
(
    const SootModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SootModel<CloudType>::~SootModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
void Foam::SootModel<CloudType>::calSoot
(
    scalarField& Y,
    const scalar& T,
    scalar& ySoot,
    scalar& nSoot,
    scalar& sootVf,
    scalar& sootIm,
    const scalar& rho
)
{}

template<class CloudType>
void Foam::SootModel<CloudType>::calSoot
(
    scalarField& Y,
    const scalar& T,
    scalar& ySoot,
    scalar& nSoot,
    scalar& sootVf,
    scalar& sootIm,
    const scalar& rho,
    scalarList& sootContribs
)
{}

template<class CloudType>
Foam::Switch Foam::SootModel<CloudType>::turbMixing() const
{
    NotImplemented;
    Info << "*** CALLING WRONG FUNCION "<<endl;
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SootModelNew.C"

// ************************************************************************* //
