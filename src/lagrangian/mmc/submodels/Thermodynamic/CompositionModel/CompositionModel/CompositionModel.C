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

#include "CompositionModel.H"
#include "reactingMixture.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel(CloudType& owner)
:
    SubModelBase<CloudType>(owner),
    thermo_(owner.thermo()),
    SLGThermo_(owner.slgThermo()),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
    (owner.thermo()).speciesData()
    ),
    mixture_("mixture", specieThermo_[0])
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type),
    thermo_(owner.thermo()),
    SLGThermo_(owner.slgThermo()),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
    (owner.thermo()).speciesData()
    ),
    mixture_("mixture", specieThermo_[0])
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const CompositionModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm),
    thermo_(cm.thermo_),
    SLGThermo_(cm.SLGThermo_),
    specieThermo_(cm.specieThermo_),
    mixture_(cm.mixture_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::~CompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::SLGThermo& Foam::CompositionModel<CloudType>::slgThermo() const
{
    return SLGThermo_;
}


template<class CloudType>
const typename CloudType::ReactionThermoType& Foam::CompositionModel<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
const Foam::basicMultiComponentMixture&
Foam::CompositionModel<CloudType>::carrier() const
{
    return slgThermo().carrier();
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::componentNames() const
{
    return slgThermo().carrier().species();
}


template<class CloudType>
const Foam::gasHThermoPhysics&
Foam::CompositionModel<CloudType>::particleMixture
(
    const scalarField& Yp
) const
{
    mixture_ = Yp[0]*specieThermo_[0];

    for (label n=1; n<specieThermo_.size(); n++)
    {
        mixture_ += Yp[n]*specieThermo_[n];
    }

    return mixture_;
}


template<class CloudType>
const Foam::scalarField Foam::CompositionModel<CloudType>::YMixture0(label patch, label patchFace) 
{
    notImplemented
    (
        "const scalarField& Foam::CompositionModel<CloudType>::YMixture0() "
        "const"
    );

    return scalarField::null();
}


template<class CloudType>
const Foam::scalarField Foam::CompositionModel<CloudType>::YMixture0(label celli) 
{
    notImplemented
    (
        "const scalarField& Foam::CompositionModel<CloudType>::YMixture0() "
        "const"
    );

    return scalarField::null();
}


template<class CloudType>
Foam::scalarField Foam::CompositionModel<CloudType>::X(const scalarField& Y) const
{
    scalarField X(Y.size());
    scalar Winv = 0.0;
    forAll(X, i)
    {
        Winv += Y[i]/specieThermo_[i].W();
        X[i] = Y[i]/specieThermo_[i].W();
    }

    tmp<scalarField> tfld = X/Winv;
    return tfld();
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::molWt(const label i) const
{
    return specieThermo_[i].W();    
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::mixtureMW(const scalarField& Y) const
{
    scalar YoverW = 0.0;
    forAll(Y, i)
    {
        YoverW += Y[i]/specieThermo_[i].W();
    }

    return 1.0 / YoverW;
}


template<class CloudType>
void Foam::CompositionModel<CloudType>::correctMassFractions()
{
    label inertIndex = -1;

    const word inertSpecie(thermo().lookup("inertSpecie"));
    
    forAllIters(this->owner(), iter)
    {
        scalar Yt = 0;
    
        forAll(iter().Y(), specieI)
        {
            if (componentNames()[specieI] != inertSpecie)
            {
                if (iter().Y()[specieI] < 0)
                    iter().Y()[specieI] = 0;
        
                Yt += iter().Y()[specieI];
            }
            else
            {
                inertIndex = specieI;
            }
        }
        
        iter().Y()[inertIndex] = scalar(1) - Yt;
        
        if (iter().Y()[inertIndex] < 0)
            iter().Y()[inertIndex] = 0;
    }
}


template<class CloudType>
const Foam::Switch
Foam::CompositionModel<CloudType>::printMoleFractionsEnabled() const
{
    return false;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CompositionModelNew.C"

// ************************************************************************* //

