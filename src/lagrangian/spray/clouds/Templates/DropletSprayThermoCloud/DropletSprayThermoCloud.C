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

#include "DropletSprayThermoCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::setModels()
{
    dropletToMMCModel_.reset
    (
        DropletToMMCModel<thermoCloudType, mmcCloudType>::New
        (
            this->subModelProperties(),
            *this
        )
    );

    liquidProperties_.reset
    (
        LiquidPropertiesModel<thermoCloudType>::New
        (
            this->subModelProperties(),
            *this
        )
    );
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::cloudReset
(
    DropletSprayThermoCloud<CloudType, mmcCloudType>& c
)
{
    CloudType::cloudReset(c);
    dropletToMMCModel_.reset(c.dropletToMMCModel_.ptr());
    liquidProperties_.reset(c.liquidProperties_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::DropletSprayThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    thermoModels_(baseParticleThermoMixture::New(thermo.thermo())),
    constProps_(this->particleProperties()),
    speciesName_(this->subModelProperties().template get<word>("dropletComponent")),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p())
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }

        const dictionary& dict = this->solution().integrationSchemes();
        word integrationType = dict.get<word>("T");

        if (integrationType == "EulerExplicit")
            integrationScheme_ = DropletSprayThermo::integrationScheme::EulerExplicit;
        else if (integrationType == "CrankNicolson")
            integrationScheme_ = DropletSprayThermo::integrationScheme::CrankNicolson;
        else
            FatalError << "No valid entry found for the integration scheme of T" << nl
                << "Valid entries are: " << nl
                << "\tEulerExplicit" << nl
                << "\tCrankNicolson" << nl 
                << exit(FatalError);
    }


    rhoTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans_" + speciesName_,
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMass, Zero)
        )
    );

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType, class mmcCloudType>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::DropletSprayThermoCloud
(
    DropletSprayThermoCloud<CloudType, mmcCloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    thermoModels_(c.thermoModels_->clone()),
    integrationScheme_(c.integrationScheme_),
    constProps_(c.constProps_),
    speciesName_(c.speciesName_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    dropletToMMCModel_(c.dropletToMMCModel_->clone()),
    liquidProperties_(c.liquidProperties_->clone())
{
    rhoTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans_" + speciesName_,
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.rhoTrans_()
        )
    );
}


template<class CloudType, class mmcCloudType>
Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::DropletSprayThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const DropletSprayThermoCloud<CloudType, mmcCloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    thermoModels_(c.thermoModels_->clone()),
    integrationScheme_(c.integrationScheme_),
    constProps_(c.constProps_),
    speciesName_(c.speciesName_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    dropletToMMCModel_(c.dropletToMMCModel_->clone()),
    liquidProperties_(c.liquidProperties_->clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
    parcel.rho() = this->liquidProperties().rhoL(this->p_[parcel.cell()],constProps_.T0());
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    parcel.mass0() = parcel.mass();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<DropletSprayThermoCloud<CloudType, mmcCloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    rhoTrans_.ref().field() = 0.0;
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::relaxSources
(
    const DropletSprayThermoCloud<CloudType, mmcCloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    typedef volScalarField::Internal dsfType;
 
    dsfType& rhoT = rhoTrans_.ref();
    const dsfType& rhoT0 = cloudOldTime.rhoTrans();
    this->relax(rhoT, rhoT0, "rho");
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::scaleSources()
{
    CloudType::scaleSources();

    typedef volScalarField::Internal dsfType;

    dsfType& rhoT = rhoTrans_.ref();
    this->scale(rhoT, "rho");
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayThermoCloud<CloudType, mmcCloudType>::info()
{
    CloudType::info();
}


// ************************************************************************* //
