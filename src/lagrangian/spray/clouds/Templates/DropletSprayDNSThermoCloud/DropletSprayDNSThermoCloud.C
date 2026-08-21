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

#include "DropletSprayDNSThermoCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::setModels()
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
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::setDistributeSource
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Create the averaging space environment
    autoPtr<aveSpace> aveVolume = aveSpace::New(dict);

    distSource_.reset
    (
        new distributeSource(mesh,aveVolume)
    );

    movingAverage_.reset
    (
        new movingAverage(mesh,distSource_().manager())
    );
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::cloudReset
(
    DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& c
)
{
    CloudType::cloudReset(c);
    liquidProperties_.reset(c.liquidProperties_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::DropletSprayDNSThermoCloud
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
    refPropertiesFactor_(this->subModelProperties().template getOrDefault<scalar>("refPropertiesFactor",0.333333)),
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
            integrationScheme_ = DropletSprayDNSThermo::integrationScheme::EulerExplicit;
        else if (integrationType == "CrankNicolson")
            integrationScheme_ = DropletSprayDNSThermo::integrationScheme::CrankNicolson;
        else
            FatalError << "No valid entry found for the integration scheme of T" << nl
                << "Valid entries are: " << nl
                << "\tEulerExplicit" << nl
                << "\tCrankNicolson" << nl 
                << exit(FatalError);

        // Read the distribution source model type
        const word distributeSourceModel = 
            this->subModelProperties().template get<word>("sourceDistribution");
        if (distributeSourceModel == "PSI")
            distSourceModelType_ = DropletSprayDNSThermo::distributeSourceModel::PSI;
        else if (distributeSourceModel == "distributed")
        {
            distSourceModelType_ = DropletSprayDNSThermo::distributeSourceModel::distributed;
            // If it is distributed, set the distribute model
            setDistributeSource
            (
                this->mesh(),
                this->subModelProperties().subDict("distributedCoeffs")
            );
        }
        else
        {
            Info << "Selecting default distribute source model: PSI" << endl;
            // Default is PSI
            distSourceModelType_ = DropletSprayDNSThermo::distributeSourceModel::PSI;
        }
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

    hsTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans_",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                // IOobject::AUTO_WRITE
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, Zero)
        )
    );

    hsCoeff_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff_",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                // IOobject::AUTO_WRITE
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, Zero)
        )
    );

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType, class mmcCloudType>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::DropletSprayDNSThermoCloud
(
    DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& c,
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
    refPropertiesFactor_(c.refPropertiesFactor_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
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

    hsTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsTrans_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsTrans_()
        )
    );
    hsCoeff_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":hsCoeff_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsCoeff_()
        )
    );
}


template<class CloudType, class mmcCloudType>
Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::DropletSprayDNSThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    thermoModels_(c.thermoModels_->clone()),
    integrationScheme_(c.integrationScheme_),
    constProps_(c.constProps_),
    speciesName_(c.speciesName_),
    refPropertiesFactor_(c.refPropertiesFactor_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    liquidProperties_(c.liquidProperties_->clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
    parcel.rho() = this->liquidProperties().rhoL(this->p_[parcel.cell()],constProps_.T0());
    parcel.Tvap() = constProps_.Tvap0();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::checkParcelProperties
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
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<DropletSprayDNSThermoCloud<CloudType,mmcCloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    rhoTrans_.ref().field() = 0.0;
    hsTrans_.ref().field() = 0.0;
    hsCoeff_.ref().field() = 0.0;
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::relaxSources
(
    const DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    typedef volScalarField::Internal dsfType;
 
    dsfType& rhoT = rhoTrans_.ref();
    const dsfType& rhoT0 = cloudOldTime.rhoTrans();
    this->relax(rhoT, rhoT0, "rho");

    dsfType& hsT = hsTrans_.ref();
    const dsfType& hsT0 = cloudOldTime.hsTrans();
    this->relax(hsT, hsT0, "hEqvE");

    dsfType& hsC = hsCoeff_.ref();
    const dsfType& hsC0 = cloudOldTime.hsCoeff();
    this->relax(hsC, hsC0, "hEqvE");
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::scaleSources()
{
    CloudType::scaleSources();

    typedef volScalarField::Internal dsfType;

    dsfType& rhoT = rhoTrans_.ref();
    this->scale(rhoT, "rho");

    dsfType& hsT = hsTrans_.ref();
    this->scale(hsT, "hEqvE");

    dsfType& hsC = hsCoeff_.ref();
    this->scale(hsC, "hEqvE");
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::distributeSources()
{
    if (distSourceModelType_ == DropletSprayDNSThermo::distributeSourceModel::distributed)
    {
        hsTrans_.ref() = distSource_->distribute(hsTrans_());
        rhoTrans_.ref() = distSource_->distribute(rhoTrans_());
        this->UTrans_.ref() = distSource_->distribute(this->UTrans_());
    }
    // If PSI do nothing
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType, class mmcCloudType>
void Foam::DropletSprayDNSThermoCloud<CloudType,mmcCloudType>::info()
{
    CloudType::info();
}


// ************************************************************************* //
