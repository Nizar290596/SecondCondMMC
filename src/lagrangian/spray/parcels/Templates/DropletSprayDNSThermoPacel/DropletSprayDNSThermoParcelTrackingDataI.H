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
#include "volFields.H" 
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Function * * * * * * * * * * * * * 

template<class ParcelType>
template<class TrackCloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::averagedField
(
    const TrackCloudType& cloud,
    const volScalarField& field
) const
{
    if (cloud.distributionModel() == DropletSprayDNSThermo::distributeSourceModel::distributed)
    {
        auto& ave = cloud.averageFunction();
        return tmp<volScalarField>
        (
            new volScalarField
            (
                field.name(),
                ave.average(field)
            )
        );
    }
    
    return field;
}



template<class ParcelType>
template<class TrackCloudType>
inline Foam::PtrList<volScalarField>
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::averageSpeciesField
(
    const TrackCloudType& cloud,
    const PtrList<volScalarField>& YField
) const
{
    if (cloud.distributionModel() == DropletSprayDNSThermo::distributeSourceModel::distributed)
    {
        PtrList<volScalarField> YFieldAveraged;
        forAll(YField,i)
        {
            // Create a new volScalarField
            YFieldAveraged.append
            (
                new volScalarField
                (
                    YField[i].name(),
                    averagedField(cloud,YField[i])
                )
            );
        }
        return YFieldAveraged;
    }
    
    return YField;
}



// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * 

template<class ParcelType>
template<class TrackCloudType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    Cp_(averagedField(cloud,cloud.thermo().thermo().Cp())),
    kappa_(averagedField(cloud,cloud.thermo().thermo().kappa())),
    Y_(averageSpeciesField(cloud,cloud.thermo().carrier().Y())),
    T_(averagedField(cloud,cloud.T())),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            T_
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
    pInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.p()
        )
    ), 
    GInterp_(nullptr),
    Tc_(Zero),
    Cpc_(Zero),
    pc_(Zero)
{   
    // const fvMesh& mesh = cloud.mesh(); 
    const speciesTable species = cloud.thermo().carrier().species();

    forAll(species,i)
    {   
        YiInterp_.emplace
        (
            species[i],
            interpolation<scalar>::New
            (
                cloud.solution().interpolationSchemes(),
                Y_[i]
            )
        );
    }

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
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::kappa() const
{
    return kappa_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::TInterp() const
{
    return *TInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::XiRInterp
(
    const word name
) const
{
    return XiRInterp_[name]();
}

template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::YiInterp
(
    const word name
) const
{
    return YiInterp_[name]();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::CpInterp() const
{
    return *CpInterp_;
}

template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::pInterp() const
{
    return *pInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::kappaInterp() const
{
    return *kappaInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::GInterp() const
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
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::Tc()
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::Cpc() const
{
    return Cpc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::Cpc()
{
    return Cpc_;
}

template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::pc() const
{
    return pc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::trackingData::pc()
{
    return pc_;
}


// ************************************************************************* //
