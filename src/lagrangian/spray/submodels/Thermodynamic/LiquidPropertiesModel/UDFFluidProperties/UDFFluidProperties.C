/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UDFFluidProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UDFFluidProperties<CloudType>::UDFFluidProperties
(
    const dictionary& dict,
    CloudType& owner
)
:
    LiquidPropertiesModel<CloudType>(dict, owner, typeName),
    saturationPressure_(Function1<scalar>::New
    (
        "saturationPressure",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    saturationTemperature_(Function1<scalar>::New
    (
        "saturationTemperature",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    latentHeat_(Function1<scalar>::New
    (
        "latentHeat",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    heatCapacityLiquid_(Function1<scalar>::New
    (
        "liquidHeatCapacity",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    absoluteEnthalpy_(Function1<scalar>::New
    (
        "liquidAbsoluteEnthalpy",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    rhoL_(Function1<scalar>::New
    (
        "rhoL",
        this->coeffDict(),
        &(owner.pMesh().thisDb()))
    ),
    binaryDiffusionCoefficient_
    (
        this->coeffDict().template get<scalar>("binaryDiffusionCoefficient")
    )
{}


template<class CloudType>
Foam::UDFFluidProperties<CloudType>::UDFFluidProperties
(
    const UDFFluidProperties<CloudType>& copyModel
)
:
    LiquidPropertiesModel<CloudType>(copyModel.owner_),
    saturationPressure_(copyModel.saturationPressure_->clone().ptr()),
    saturationTemperature_(copyModel.saturationTemperature_->clone().ptr()),
    latentHeat_(copyModel.latentHeat_->clone().ptr()),
    heatCapacityLiquid_(copyModel.heatCapacityLiquid_->clone().ptr()),
    absoluteEnthalpy_(copyModel.absoluteEnthalpy_->clone().ptr()),
    rhoL_(copyModel.rhoL_->clone().ptr()),
    binaryDiffusionCoefficient_(copyModel.binaryDiffusionCoefficient_)
{}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::pSat
(
    const scalar T
) const
{
    return saturationPressure_->value(T);
}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::TSat
(
    const scalar p
) const
{
    return saturationTemperature_->value(p);
}



template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::Df
(
    const scalar p, 
    const scalar T,
    const scalar Wb
) const
{
    // Wb is here ignored
    return binaryDiffusionCoefficient_*std::pow(T,1.75)/p;
}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::Lv
(
    const scalar p,
    const scalar T
) const
{
    return latentHeat_->value(T);
}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::CpL
(
    const scalar p, 
    const scalar T
) const
{
    return heatCapacityLiquid_->value(T);
}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::rhoL
(
    const scalar p, 
    const scalar T
) const
{
    return rhoL_->value(T);
}


template<class CloudType>
Foam::scalar Foam::UDFFluidProperties<CloudType>::Ha
(
    const scalar p,
    const scalar T
) const
{
    return absoluteEnthalpy_->value(T);
}

// ************************************************************************* //
