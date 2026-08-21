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
inline const Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType, class mmcCloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::constProps() const
{
    return constProps_;
}
 
 
template<class CloudType, class mmcCloudType>
inline typename CloudType::particleType::constantProperties&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::constProps()
{
    return constProps_;
}
 

template<class CloudType, class mmcCloudType>
inline scalar Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::refPropertiesFactor() const
{
    return refPropertiesFactor_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::SLGThermo& Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::thermo() const
{
    return thermo_;
}

template<class CloudType, class mmcCloudType>
inline const Foam::baseParticleThermoMixture& 
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::thermoModels() const
{
    return thermoModels_();
}


template<class CloudType, class mmcCloudType>
inline const Foam::LiquidPropertiesModel
<
    Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>
>& 
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::liquidProperties() const
{
    return liquidProperties_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::volScalarField& Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::T() const
{
    return T_;
}
 
 
template<class CloudType, class mmcCloudType>
inline const Foam::volScalarField& Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::p() const
{
    return p_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::word& 
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::dropletSpeciesName() const
{
    return speciesName_;
}


template<class CloudType, class mmcCloudType>
inline const Foam::movingAverage&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::averageFunction() const
{
    return movingAverage_();
}


template<class CloudType, class mmcCloudType>
inline const Foam::DropletSprayDNSThermo::distributeSourceModel&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::distributionModel() const
{
    return distSourceModelType_;
}


template<class CloudType, class mmcCloudType>
void
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::setMMCCloud
(
    mmcCloudType& cloud
)
{
    mmcCloud_.ref(cloud);
}


template<class CloudType, class mmcCloudType>
inline const mmcCloudType&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::mmcCloud() const
{
    return mmcCloud_.ref();
}


template<class CloudType, class mmcCloudType>
inline mmcCloudType&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::mmcCloud()
{
    return mmcCloud_.ref();
}


template<class CloudType, class mmcCloudType>
inline const DropletToMMCModel<DropletSprayDNSThermoCloud<CloudType,mmcCloudType>, mmcCloudType>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::dropletToMMCCoupling() const
{
    return dropletToMMCModel_();
}

template<class CloudType, class mmcCloudType>
inline const Foam::DropletSprayDNSThermo::integrationScheme&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::integrationScheme() const
{
    return integrationScheme_;
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::Ep() const
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
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::ap() const
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
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::sigmap() const
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
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::rhoTrans() const
{
    return rhoTrans_();
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::rhoTrans()
{
    return rhoTrans_.ref();
}


template<class CloudType, class mmcCloudType>
inline
const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::hsTrans() const
{
    return hsTrans_();
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::hsTrans()
{
    return hsTrans_.ref();
}


template<class CloudType, class mmcCloudType>
inline
const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::hsCoeff() const
{
    return hsCoeff_();
}
 
 
template<class CloudType, class mmcCloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::hsCoeff()
{
    return hsCoeff_.ref();
}
 

template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::Srho() const
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

// Sources for DNS filtered field
// ==============================

template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::SrhoDNS
(
    const fvMesh& DNSMesh,
    const filterDNSField& filter
) const
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
            DNSMesh,
            dimensionedScalar
            (
                rhoTrans_().dimensions()/dimTime/dimVolume, Zero
            )
        )
    );
 
    if (this->solution().coupled())
    {
        scalarField& sourceField = trhoTrans.ref();

        if (sourceField.size() != rhoTrans_().size())
        {
            tmp<volScalarField::Internal> rhoTransDNS = 
                filter.mapSourceFilterFieldToDNS(rhoTrans_());
            sourceField += rhoTransDNS;
        }
        else
        {
            sourceField += rhoTrans_();
        }

        sourceField /= this->db().time().deltaTValue()*DNSMesh.V();
    }
 
    return trhoTrans;
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::SrhoDNS
(
    const fvMesh& DNSMesh,
    const filterDNSField& filter,
    volScalarField& rho
) const
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
                DNSMesh,
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
            )
        );
 
        scalarField& sourceField = trhoTrans.ref();
 
        if (this->solution().semiImplicit("rho"))
        {
            if (rho.size() != rhoTrans().size())
            {
                tmp<volScalarField::Internal> rhoTransDNS = 
                    filter.mapSourceFilterFieldToDNS(rhoTrans());
                sourceField += rhoTransDNS();
                sourceField /= this->db().time().deltaTValue()*DNSMesh.V();
            }
            else
            {
                sourceField += rhoTrans_();
                sourceField /= this->db().time().deltaTValue()*DNSMesh.V();
            }
            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();
            if (rho.size() != rhoTrans().size())
            {
                tmp<volScalarField::Internal> rhoTransDNS = 
                    filter.mapSourceFilterFieldToDNS(rhoTrans());
                sourceField += rhoTransDNS();
            }
            else
            {
                sourceField += rhoTrans_();
            }
 
            fvm.source() = -trhoTrans()/this->db().time().deltaT();
 
            return tfvm;
        }
    }
 
    return tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::SYFuDNS
(
    const fvMesh& DNSMesh,
    const filterDNSField& filter,
    volScalarField& Yi
) const
{
    
    if (this->solution().coupled() && Yi.name() == "C2H5OH")
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
                DNSMesh,
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
            )
        );
 
        scalarField& sourceField = trhoTrans.ref();
 
        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        if (sourceField.size() != rhoTrans_().size())
        {
            tmp<volScalarField::Internal> rhoTransDNS = 
                filter.mapSourceFilterFieldToDNS(rhoTrans_());
            sourceField += rhoTransDNS;
        }
        else
        {
            sourceField += rhoTrans_();
        }

        fvm.source() = -trhoTrans()/this->db().time().deltaT();

        return tfvm;
    }

    return tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::ShDNS
(
    const fvMesh& DNSMesh,
    const filterDNSField& filter,
    volScalarField& hs
) const
{
    if (debug)
    {
        Info<< "hsTrans min/max = " << min(hsTrans()).value() << ", "
            << max(hsTrans_()).value() << nl
            << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
            << max(hsCoeff_()).value() << endl;
    }

    if (this->solution().coupled())
    {

        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        if (hs.size() != this->hsTrans().size())
        {
            // Map the hsTrans field to the DNS mesh
            tmp<volScalarField::Internal> hsTransDNS = 
                filter.mapSourceFilterFieldToDNS(hsTrans_());
            const volScalarField::Internal& refer = hsTransDNS();
            fvm.source() = -refer/(this->db().time().deltaT());

             Info<< "hsTrans min/max = " << min(hsTrans_()).value() << ", "
            << max(hsTrans_()).value() << endl;
            Info<< "hsTransDNS min/max = " << min(hsTransDNS()).value() << ", "
            << max(hsTransDNS()).value() << endl;
        }
        else
        {
            Info<< "hsTrans min/max = " << min(hsTrans_()).value() << ", "
            << max(hsTrans_()).value() << endl;
            fvm.source() = -hsTrans_()/(this->db().time().deltaT());
        }
        return tfvm;
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(hs, dimEnergy/dimTime));
}



template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::SUDNS
(
    const fvMesh& DNSMesh,
    const filterDNSField& filter,
    volVectorField& U
) const
{
    if (debug)
    {
        Pout<< "UTrans min/max = " << min(this->UTrans()).value() << ", "
            << max(this->UTrans()).value() << nl
            << "UCoeff min/max = " << min(this->UCoeff()).value() << ", "
            << max(this->UCoeff()).value() << endl;
    }
 
    dimensionSet dim(dimForce);
 
    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("U"))
        {
            volScalarField::Internal
                Vdt(DNSMesh.V()*this->db().time().deltaT());
 
            if (U.size() != this->UTrans().size())
            {
                tmp<volVectorField::Internal> UTransDNS = 
                    filter.mapSourceFilterFieldToDNS(this->UTrans());
                tmp<volScalarField::Internal> UCoeffDNS = 
                    filter.mapSourceFilterFieldToDNS(this->UCoeff());
                return UTransDNS()/Vdt - fvm::Sp(UCoeffDNS()/Vdt, U) + UCoeffDNS()/Vdt*U;
            }
            else
            {
                return this->UTrans()/Vdt - fvm::Sp(this->UCoeff()/Vdt, U) + this->UCoeff()/Vdt*U;
            }
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dim));
            fvVectorMatrix& fvm = tfvm.ref();
            
            if (U.size() != this->UTrans().size())
            {
                tmp<volVectorField::Internal> UTransDNS = 
                    filter.mapSourceFilterFieldToDNS(this->UTrans());
                fvm.source() = -UTransDNS/(this->db().time().deltaT());
            }
            else
            {
                fvm.source() = -this->UTrans()/(this->db().time().deltaT());
            }

            return tfvm;
        }
    }
 
    return tmp<fvVectorMatrix>::New(U, dim);
}

// --- END of DNS Sources 

template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::Sh(volScalarField& hs) const
{
    if (debug)
    {
        Info<< "hsTrans min/max = " << min(hsTrans()).value() << ", "
            << max(hsTrans_()).value() << nl
            << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
            << max(hsCoeff_()).value() << endl;
    }

    if (this->solution().coupled())
    {

        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        fvm.source() = -hsTrans_()/(this->db().time().deltaT());

        return tfvm;

    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(hs, dimEnergy/dimTime));
}


template<class CloudType, class mmcCloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::Srho(volScalarField& rho) const
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
