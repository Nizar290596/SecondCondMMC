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

#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
inline const Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType, class mmcCloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::constProps() const
{
    return constProps_;
}
 
 
template<class CloudType, class mmcCloudType>
inline typename CloudType::particleType::constantProperties&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::constProps()
{
    return constProps_;
}
 
 
template<class CloudType, class mmcCloudType>
inline const Foam::SLGThermo& Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::thermo() const
{
    return thermo_;
}

template<class CloudType, class mmcCloudType>
inline const Foam::baseParticleThermoMixture& 
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::thermoModels() const
{
    return thermoModels_();
}


template<class CloudType, class mmcCloudType>
inline const Foam::LiquidPropertiesModel
<
    Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>
>& 
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::liquidProperties() const
{
    return liquidProperties_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::volScalarField& Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::T() const
{
    return T_;
}
 
 
template<class CloudType, class mmcCloudType>
inline const Foam::volScalarField& Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::p() const
{
    return p_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::word& 
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::dropletSpeciesName() const
{
    return speciesName_;
}


template<class CloudType, class mmcCloudType>
void
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::setMMCCloud
(
    mmcCloudType& cloud
)
{
    mmcCloud_.ref(cloud);
}


template<class CloudType, class mmcCloudType>
inline const mmcCloudType&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::mmcCloud() const
{
    return mmcCloud_.ref();
}


template<class CloudType, class mmcCloudType>
inline mmcCloudType&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::mmcCloud()
{
    return mmcCloud_.ref();
}


template<class CloudType, class mmcCloudType>
inline const DropletToMMCModel<DropletSprayThermoCloud<CloudType,mmcCloudType>, mmcCloudType>&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::dropletToMMCCoupling() const
{
    return dropletToMMCModel_();
}


template<class CloudType, class mmcCloudType>
inline const Foam::DropletSprayThermo::integrationScheme&
Foam::DropletSprayThermoCloud<CloudType,mmcCloudType>::integrationScheme() const
{
    return integrationScheme_;
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );
 
 
    return tEp;
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );
 
 
    return tap;
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );
 
 
    return tsigmap;
}


template<class CloudType, class mmcCloudType>
inline
const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::rhoTrans() const
{
    return rhoTrans_();
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::rhoTrans()
{
    return rhoTrans_.ref();
}

 
template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::Srho() const
{
    tmp<volScalarField::Internal> trhoTrans
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_().dimensions()/dimTime/dimVolume, Zero
            )
        )
    );
 
    if (this->solution().coupled())
    {
        scalarField& sourceField = trhoTrans.ref();
        sourceField += rhoTrans_();
        sourceField /= this->db().time().deltaTValue()*this->mesh().V();
    }
 
    return trhoTrans;
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        tmp<volScalarField> trhoTrans
        (
            new volScalarField
            (
                IOobject
                (
                    this->name() + ":rhoTrans",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
            )
        );
 
        scalarField& sourceField = trhoTrans.ref();
 
        if (this->solution().semiImplicit("rho"))
        {
            sourceField += rhoTrans_();
            sourceField /= this->db().time().deltaTValue()*this->mesh().V();
 
            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();
 
                sourceField += rhoTrans_();
 
            fvm.source() = -trhoTrans()/this->db().time().deltaT();
 
            return tfvm;
        }
    }
 
    return tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
}


// ************************************************************************* //
