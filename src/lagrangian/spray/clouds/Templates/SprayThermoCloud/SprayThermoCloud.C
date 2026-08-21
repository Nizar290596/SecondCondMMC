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

#include "SprayThermoCloud.H"
#include "ThermoParcel.H"

#include "EvaporationModel.H"
#include "EnvelopeModel.H"
#include "CpModel.H"
#include "FPRadiationModel.H"
#include "basicAerosolReactingPopeCloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::setModels()
{
    evaporationModel_.reset
    (
        EvaporationModel<SprayThermoCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    envelopeModel_.reset
    (
        EnvelopeModel<SprayThermoCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    cpModel_.reset
    (
        CpModel<SprayThermoCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    fpRadiationModel_.reset
    (
        FPRadiationModel<SprayThermoCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "TD",
            this->solution().integrationSchemes()
        ).ptr()
    );

    MIntegrator_.reset
    (
        integrationScheme::New
        (
            "mD",
            this->solution().integrationSchemes()
        ).ptr()
    );
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::cloudReset(SprayThermoCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    evaporationModel_.reset(c.evaporationModel_.ptr());
    envelopeModel_.reset(c.envelopeModel_.ptr());
    cpModel_.reset(c.cpModel_.ptr());    
    fpRadiationModel_.reset(c.fpRadiationModel_.ptr());    
    TIntegrator_.reset(c.TIntegrator_.ptr());
    MIntegrator_.reset(c.MIntegrator_.ptr());
}

template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::parcelIdAddOne(const label id)
{
    label add = 1;
    if(this->parcelId().size() < 1.0)
    {
        this->parcelId().append(id);
    }
    else
    {
        forAll(parcelId(), index)
        {
            if(parcelId()[index] == id)
            {
                add = 0;            
            }
        }
        if(add == 1)
        {
            this->parcelId().append(id);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayThermoCloud<CloudType>::SprayThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& muEff,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    const volScalarField& f,
    basicAerosolReactingPopeCloud& pope,
    const volScalarField& magSqrGradf,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        muEff,
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties(), this->solution().active()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    evaporationModel_(nullptr),
    envelopeModel_(nullptr), 
    cpModel_(nullptr),
    fpRadiationModel_(nullptr),
    TIntegrator_(nullptr),
    MIntegrator_(nullptr),
    f_(f),
    magSqrGradf_(magSqrGradf),
    rand_f_Gen_(Pstream::myProcNo()+1),
    pope_(pope),
    specieId_(),
    nFuelSpecies_(1),
    linkFG_(3),
    nProductSpecies_(1),
    linkPG_(3),
    YProducts_(3),
    MTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + "MTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedVector("zero", dimMass*dimVelocity, Zero)
        )
    ),
    MTransRate_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "MTransRate",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
        )
    ),
    emptySuperCell_(0)
{
    if (this->solution().active())
    {
        setModels();
        setEulerianStatistics();
        setFuelMix();
        setProductMix();
        if (readFields)
        {
            if(this->size()>0)
            {
                parcelType::readFields(*this, pope_.composition()); 
            }
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::SprayThermoCloud<CloudType>::SprayThermoCloud
(
    SprayThermoCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    evaporationModel_(c.evaporationModel_->clone()),
    envelopeModel_(c.envelopeModel_->clone()),
    cpModel_(c.cpModel_->clone()),
    fpRadiationModel_(c.fpRadiationModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    MIntegrator_(c.MIntegrator_->clone()),
    f_(c.f()),
    magSqrGradf_(c.magSqrGradf()),
    rand_f_Gen_(c.rand_f_Gen_),
    pope_(c.pope()),
    specieId_(c.specieId()),
    nFuelSpecies_(c.nFuelSpecies()),
    linkFG_(c.linkFG()),
    nProductSpecies_(c.nProductSpecies()),
    linkPG_(c.linkPG()),
    YProducts_(c.YProducts()),
    MTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + "MTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.MTrans()
        )
    ),
    MTransRate_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "MTransRate",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.MTransRate()
        )
    )
{}


template<class CloudType>
Foam::SprayThermoCloud<CloudType>::SprayThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const SprayThermoCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    evaporationModel_(nullptr),
    envelopeModel_(nullptr),
    cpModel_(nullptr),    
    fpRadiationModel_(nullptr),    
    TIntegrator_(nullptr),
    MIntegrator_(nullptr),
    f_(c.f()),
    magSqrGradf_(c.magSqrGradf()),
    rand_f_Gen_(c.rand_f_Gen_),
    pope_(c.pope()),
    specieId_(c.specieId()),
    nFuelSpecies_(c.nFuelSpecies()),
    linkFG_(c.linkFG()),
    nProductSpecies_(c.nProductSpecies()),
    linkPG_(c.linkPG()),
    YProducts_(c.YProducts()),
    MTrans_(nullptr),
    MTransRate_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayThermoCloud<CloudType>::~SprayThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//link the gas thermo species id to the fuel components
template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::setFuelMix()
{
    nFuelSpecies_ = this->particleProperties().subDict("composition").size();
    linkFG_.setSize(nFuelSpecies(), 0.0);

    label nfs = 0;
    forAll(pope_.composition().componentNames(), ns)
    {
        word fuelName = pope_.composition().componentNames()[ns];
        if(this->particleProperties().subDict("composition").found(fuelName))
        {
            linkFG_[nfs] = ns;
            nfs++;
        }
    }
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::setProductMix()
{
    nProductSpecies_ = this->particleProperties().subDict("productComposition").size();
    linkPG_.setSize(nProductSpecies(), 0.0);
    YProducts_.setSize(nProductSpecies(), 0.0);

    label nfs = 0;
    forAll(pope_.composition().componentNames(), ns)
    {
        word productName = pope_.composition().componentNames()[ns];
        if(this->particleProperties().subDict("productComposition").found(productName))
        {
            linkPG_[nfs] = ns;
            YProducts_[nfs] = readScalar(this->particleProperties().subDict("productComposition").lookup(productName));
Info << "Product: " << productName << ", value: " << YProducts_[nfs] << endl;
            nfs++;
        }
    }
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.T() = constProps_.T0();
    parcel.Cp() = cpModel().Cp(constProps_.T0()); 
    parcel.nThermoSpecies() = pope_.composition().componentNames().size();
    parcel.molWtSstate() = constProps_.molWtSstate0();
    parcel.nFuelSpecies() = this->nFuelSpecies();
    parcel.YFuel().setSize(this->nFuelSpecies(), 0.0);

    label nfs = 0;
    forAll(pope_.composition().componentNames(), ns)
    {
        word fuelName = pope_.composition().componentNames()[ns];
        if(this->particleProperties().subDict("composition").found(fuelName))
        {
            parcel.YFuel()[nfs] = readScalar(this->particleProperties().subDict("composition").lookup(fuelName));
	    nfs++;
	}
    }    
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<SprayThermoCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    MTrans_->field() = Zero;
    MTransRate_->field() = 0.0;
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::relaxSources
(
    const SprayThermoCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(MTrans_(), cloudOldTime.MTrans(), "U");
    this->relax(MTransRate_(), cloudOldTime.MTransRate(), "rho");
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(MTrans_(), "U");
    this->scale(MTransRate_(), "rho");
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);
        
        this->solve(*this,td);
    }

    emptySuperCell().clear();
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::setEulerianStatistics()
{
    if (this->eulerianStatsDict().found("mV"))
    {
        const dimensionSet dim = dimMass;        
        this->eulerianStats().newProperty("mV",dim);
    }    

    if (this->eulerianStatsDict().found("fGas"))
    {
        const dimensionSet dim = dimless;        
        this->eulerianStats().newProperty("fGas",dim);
    }  

    if (this->eulerianStatsDict().found("TD"))
    {
        const dimensionSet dim = dimTemperature;        
        this->eulerianStats().newProperty("TD",dim);
    }   

    if (this->eulerianStatsDict().found("TGas"))
    {
        const dimensionSet dim = dimTemperature;        
        this->eulerianStats().newProperty("TGas",dim);
    }

    if (this->eulerianStatsDict().found("mRate"))
    {  
        this->eulerianStats().newProperty("mRate",dimensionSet(1, 0, -1, 0, 0));
    }
}

template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();       
           
    forAllIters(*this, iter)
    {   
        this->eulerianStats().findCell(iter().position());  
        
        if (this->eulerianStatsDict().found("fGas"))  
        {      
            scalar wt = 1.0;                 
            this->eulerianStats().calculate("fGas", wt, iter().fGas()); 
        }

        if (this->eulerianStatsDict().found("TD"))  
        {      
            scalar wt = 1.0;                 
            this->eulerianStats().calculate("TD", wt, iter().T()); 
        }

        if (this->eulerianStatsDict().found("TGas"))  
        {      
            scalar wt = 1.0;                 
            this->eulerianStats().calculate("TGas", wt, iter().TGas()); 
        }

        if (this->eulerianStatsDict().found("mRate"))  
        {      
            scalar wt = 1.0;                 
            this->eulerianStats().calculate("mRate", wt, iter().mFlux() * iter().d() * iter().d() * 3.14); 
        }
    }
}

/*
template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::gatherHeader(DynamicList<word> &header)
{
    CloudType::gatherHeader(header);
      
    header.append("fGas");

    header.append("fSurf");
        
    header.append("TD");

    header.append("TGas");
}


template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::gatherParcelData
(
    List<DynamicList<scalar> > &parcelSampleData
)
{
    CloudType::gatherParcelData(parcelSampleData);
    
    label nfGas    = 0;
    label nfSurf = nfGas + 1;
    label nTD = nfSurf + 1;
    label nTGas = nTD + 1;
    label numP   = nTGas + 1;
    
    label i = 0;
    //- list length is number of physical properties plus number of PSD properties 
    forAllIters(*this, iter)
    {
        DynamicList<scalar> p(numP, 0.0);

        p[nfGas] = iter().fGas();

        p[nfSurf] = iter().fSurf();

        p[nTD] = iter().T();

        p[nTGas] = iter().TGas();
                
        parcelSampleData[i].append(p);
        i++;
    }   
}
*/
template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::info()
{
    CloudType::info();

    Info<< "    Temperature min/max             = " << Tmin() << ", " << Tmax()
        << endl;
    Info<< "    No. of empty super cell(s)      = " << emptySuperCell().size() << endl;

    scalar d32 = 1.0e+6*this->Dij(3, 2);
    scalar d10 = 1.0e+6*this->Dij(1, 0);
    scalar dMax = 1.0e+6*this->Dmax();

    Info << "    D10, D32, Dmax (mu)             = " << d10 << ", " << d32
         << ", " << dMax << nl;
}

template<class CloudType>
void Foam::SprayThermoCloud<CloudType>::writeFields() const
{
    if (this->size())
    {
        CloudType::parcelType::writeFields(*this, pope_.composition());
    }
}
// ************************************************************************* //
