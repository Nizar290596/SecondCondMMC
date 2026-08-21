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

#include "OFLiquidProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::OFLiquidProperties<CloudType>::OFLiquidProperties
(
    const dictionary& dict,
    CloudType& owner
)
:
    LiquidPropertiesModel<CloudType>(dict, owner, typeName),
    liqProp_
    (
        liquidProperties::New
        (
            this->coeffDict().template get<word>("fluidName")
        )
    )
{}


template<class CloudType>
Foam::OFLiquidProperties<CloudType>::OFLiquidProperties
(
    const OFLiquidProperties<CloudType>& copyModel
)
:
    LiquidPropertiesModel<CloudType>(copyModel.owner_),
    liqProp_(copyModel.liqProp_->clone().ptr())
{}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::pSat
(
    const scalar T
) const
{
    // Note: for the saturation pressure the pressure is irrelevant
    //       use dummy value here
    return liqProp_->pv(1E+5,T);
}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::TSat
(
    const scalar p
) const
{
    return liqProp_->pvInvert(p);
}



template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::Df
(
    const scalar p, 
    const scalar T,
    const scalar Wb
) const
{
    return liqProp_->D(p,T,Wb);
}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::Lv
(
    const scalar p, 
    const scalar T
) const
{
    return liqProp_->hl(p,T);
}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::CpL
(
    const scalar p, 
    const scalar T
) const
{
    return liqProp_->Cp(p,T);
}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::rhoL
(
    const scalar p, 
    const scalar T
) const
{
    return liqProp_->rho(p,T);
}


template<class CloudType>
Foam::scalar Foam::OFLiquidProperties<CloudType>::Ha
(
    const scalar p,
    const scalar T
) const
{
    return liqProp_->h(p,T);
}

// ************************************************************************* //
