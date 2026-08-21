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

#include "mmcVarSet.H"

template<class ParcelType>
template<class TrackCloudType>
inline Foam::DropletSprayThermoParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    Cp_(cloud.thermo().thermo().Cp()),
    kappa_(cloud.thermo().thermo().kappa()),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    CpInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            Cp_
        )
    ),
    kappaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            kappa_
        )
    ),
    GInterp_(nullptr),
    Tc_(Zero),
    Cpc_(Zero)
{
    const auto& mmcCloud = cloud.mmcCloud();
    if (mmcCloud.mixing().numXiR())
    {
        const mmcVarSet& setOfXi(mmcCloud.mixing().XiR());        
        const labelHashTable& XiIndexes = setOfXi.rVarInXi();

        for (const word& XiRName : mmcCloud.mixing().XiRNames())
        {
            if (setOfXi.Vars(XiIndexes[XiRName]).refType()=="interpolated")
            {                   
                XiRInterp_.emplace
                (
                    XiRName,
                    interpolation<scalar>::New
                    (
                        mmcCloud.solution().interpolationSchemes(),
                        setOfXi.Vars(XiIndexes[XiRName]).field()
                    )
                );
            }
        }
    }
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::kappa() const
{
    return kappa_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::TInterp() const
{
    return *TInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::XiRInterp
(
    const word name
) const
{
    return XiRInterp_[name]();
}



template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::CpInterp() const
{
    return *CpInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::kappaInterp() const
{
    return *kappaInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayThermoParcel<ParcelType>::trackingData::GInterp() const
{
    if (!GInterp_)
    {
        FatalErrorInFunction
            << "Radiation G interpolation object not set"
            << abort(FatalError);
    }

    return *GInterp_;
}


template<class ParcelType>
inline Foam::scalar Foam::DropletSprayThermoParcel<ParcelType>::trackingData::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayThermoParcel<ParcelType>::trackingData::Tc()
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::DropletSprayThermoParcel<ParcelType>::trackingData::Cpc() const
{
    return Cpc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayThermoParcel<ParcelType>::trackingData::Cpc()
{
    return Cpc_;
}


// ************************************************************************* //
