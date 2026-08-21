/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ItoPopeCloud

Description
    LES-FDF cloud

SourceFiles
    ItoPopeCloudI.H
    ItoPopeCloud.C
    ItoPopeCloudIO.C

Author
Edison GE
y.ge1222@gmail.com

\*---------------------------------------------------------------------------*/
#include "ItoPopeCloud.H"
#include "InflowBoundaryModel.H"
#include "interpolation.H"
#include "fvMesh.H"
#include "boundBox.H"
#include "fvc.H"
#include "turbulentFluidThermoModel.H"
#include "passiveParticleCloud.H"
#include "meshSearch.H"


namespace Foam
{
    // The following functions are not specific to class ItoPopeCloud and should 
    // be replace by their own .H and .C files. It looks like original is 
    // mapFields.

    static const scalar perturbFactor = 1E-6;

    const Foam::dimensionedScalar SMALL_WT("SMALL_WT", Foam::dimless, Foam::SMALL);
    
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ItoPopeCloud<CloudType>::ItoPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const Switch initAtCnstr,
    bool readFields
)
    : 
    CloudType
    (
        mesh,
        cloudName,
        false
    ),
    
    itoPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    mesh_(mesh), 
     
    cloudProperties_
    (
        IOobject
        (
            createCloudPropertiesName(cloudName),
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    solution_
    (
        cloudProperties_.subDict("solution")
    ),   
    
    subModelProperties_
    (
        cloudProperties_.subDict("subModels")
    ),
    
    writeMassConsistency_
    (
        cloudProperties_.subDict("particleManagement").getOrDefault
        (
            "writeMassConsistencyCheck",
            false
        )
    ),

    pManager_(mesh,cloudProperties_.subDict("particleManagement")),
    
    runTime_
    (
        U.time()
    ),
    
    time0_
    (
        runTime_.value()
    ),
    
    turbModelPtr_(nullptr),
    
    rndGen_(Pstream::myProcNo()+1),
    
    U_(U),
   
    DEff_(DEff),
    
    rho_(rho),

    gradRho_(gradRho),

    cellLengthScale_(mag(cbrt(mesh_.V()))),
    
    eulerianStatsDict_
    (
        cloudProperties_.subDict("eulerianStatistics")
    ),
    
    eulerianStats_(mesh_,eulerianStatsDict_),
        
    inflowBoundaryModel_(nullptr),

    fvMinusPpMass_(0.0)
{
    Info << "Creating Ito Pope Particle Cloud." << nl << endl;

    setModels();
    
    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Ito Pope particle cloud data from file." << endl;
        
            particleType::itoParticleIOType::readFields(*this);
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of Ito Pope particles into the finite volume field." << nl << endl;
            
                initReleaseParticles();
            }
        }
    }
    
    //-Text file for debugging
    if (writeMassConsistency() && Pstream::myProcNo() == 0)
    {
        massConsistencyCheck_.reset
        (   
            new OFstream("massConsistencyCheck.dat")
        );
        massConsistencyCheck_() << "Text file for debugging." << nl;
    }
}



template<class CloudType>
Foam::ItoPopeCloud<CloudType>::ItoPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const Switch initAtCnstr,
    bool readFields
)
    : 
    CloudType
    (
        mesh,
        cloudName,
        false
    ),
    
    itoPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    mesh_(mesh), 
     
    cloudProperties_
    (
        IOobject
        (
            createCloudPropertiesName(cloudName),
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    solution_
    (
        cloudProperties_.subDict("solution")
    ),   
    
    subModelProperties_
    (
        cloudProperties_.subDict("subModels")
    ),
    
    writeMassConsistency_
    (
        cloudProperties_.subDict("particleManagement").getOrDefault
        (
            "writeMassConsistencyCheck",
            false
        )
    ),

    pManager_(mesh,cloudProperties_.subDict("particleManagement")),
    
    runTime_
    (
        mesh.time()
    ),
    
    time0_
    (
        runTime_.value()
    ),
    
    turbModelPtr_(nullptr),
    
    rndGen_(Pstream::myProcNo()+1),
    
    U_(lookupOrConstruct<vector>("U",dimVelocity)),
   
    DEff_(lookupOrConstruct<scalar>("DEff",dimViscosity)),
    
    rho_(lookupOrConstruct<scalar>("rho",dimDensity)),

    gradRho_(lookupOrConstruct<vector>("gradRho",dimDensity*dimViscosity/dimLength)),

    cellLengthScale_(mag(cbrt(mesh_.V()))),
    
    eulerianStatsDict_
    (
        cloudProperties_.subDict("eulerianStatistics")
    ),
    
    eulerianStats_(mesh_,eulerianStatsDict_),
        
    inflowBoundaryModel_(nullptr),

    fvMinusPpMass_(0.0)
{
    Info << "Creating Ito Pope Particle Cloud." << nl << endl;

    setModels();
    
    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Ito Pope particle cloud data from file." << endl;
        
            particleType::itoParticleIOType::readFields(*this);
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of Ito Pope particles into the finite volume field." << nl << endl;
            
                initReleaseParticles();
            }
        }
    }
    
    //-Text file for debugging
    if (writeMassConsistency() && Pstream::myProcNo() == 0)
    {
        massConsistencyCheck_.reset
        (   
            new OFstream("massConsistencyCheck.dat")
        );
        massConsistencyCheck_() << "Text file for debugging." << nl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ItoPopeCloud<CloudType>::~ItoPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::word Foam::ItoPopeCloud<CloudType>::createCloudPropertiesName
(
    const word& cloudName
) const
{
    // try to split the cloudName at the seperator "_"
    std::stringstream ss(cloudName);
    std::string buffer;
    std::getline(ss,buffer,'_');
    std::string cloudNo;
    std::getline(ss,cloudNo);
    if (!cloudNo.empty())
    {
        // Convert to number
        label cloudNoLabel = std::stoi(cloudNo);
        if (cloudNoLabel > 1)
            return "cloudProperties_"+cloudNo;
        else
            return "cloudProperties";
    }
    return "cloudProperties";
}




template<class CloudType>
template<class Type>
Foam::GeometricField<Type, fvPatchField, volMesh>&
Foam::ItoPopeCloud<CloudType>::lookupOrConstruct
(
    const word fieldName, 
    dimensionSet units
)
{
    typedef Foam::GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& mesh = this->mesh_;

    if (!mesh.objectRegistry::foundObject<volFieldType>(fieldName))
    {
        WarningInFunction 
            << fieldName << " could not be found in objectRegistry "
            << "and is now constructed" << endl;

        volFieldType* fPtr
        (
            new volFieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		        dimensioned<Type>("tmp", units, pTraits<Type>::zero)
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return 
    (
        mesh.objectRegistry::lookupObjectRef<volFieldType>(fieldName)
    );
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::setModels()
{
    inflowBoundaryModel_.reset
    (
        InflowBoundaryModel<ItoPopeCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
    
    const objectRegistry& obr = this->mesh();
    const word turbName = "turbulenceProperties";//"LESProperties";//

    if (obr.foundObject<Foam::compressible::LESModel>(turbName))
    {
        turbModelPtr_= &obr.lookupObject<Foam::compressible::LESModel>(turbName);
    }
    else
    {

// Commented to make it work for RANS models
/*
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::ItoPopeCloud<CloudType>::setModels()"
        )
        << "Turbulence model not found in mesh database" << nl
        << "Database objects include: " << obr.sortedToc()
        << abort(FatalError);
*/
    }
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::cloudReset(ItoPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    rndGen_ = c.rndGen_;

    inflowBoundaryModel_.reset(c.inflowBoundaryModel_.ptr());
}

 
template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::initReleaseParticles()
{    
    List<particleNumberController::pContainer> listParticlePosMass = 
        pManager_.initParticlePosBasedOnField(rho_);
    
    //- Unpack the released particle list and create particles of the correct type
    label numInserted = 0;
    
    for (auto& p : listParticlePosMass)
    {        
        addNewParticle
        (
            p.location,
            p.cell,
            p.tetFace,
            p.tetPt,
            p.mass/deltaM_,
            p.mass
        );
        
        numInserted++;
    }
    
    reduce(numInserted, sumOp<label>());
    
    Info << "Number of Pope particles released into domain: " << numInserted << endl;
}


template<class CloudType>
template<class TrackCloudType>
void Foam::ItoPopeCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename particleType::trackingData& td
)
{
    //- Cache the value of deltaT for this timestep
    storeDeltaT();
    
    td.resetCounters();

    this->move(cloud,td);

    this->updateEulerianStatistics();
    
    label numParticles=this->size();
    td.reportCounters(numParticles);
    
    //- Particle weight & number control
    pManager_.correct(*this);

    //- Check mass consistency
    scalar fvMass = 0;
    scalar ppMass = 0;
    
    this->checkMassConsistency(fvMass,ppMass);

    fvMinusPpMass() = fvMass - ppMass;

    reduce(fvMass, sumOp<scalar>());
    reduce(ppMass, sumOp<scalar>());
    Info << "\tMass in fv domain is " << fvMass 
         << ". Mass in particle field is "<< ppMass << endl;
         
    if (writeMassConsistency() && Pstream::myProcNo() == 0)
        massConsistencyCheck_() << fvMass << token::TAB << ppMass << token::TAB << numParticles << endl;
}


template<class CloudType>
template<class TrackCloudType>
void Foam::ItoPopeCloud<CloudType>::move
(
    TrackCloudType& cloud,
    typename particleType::trackingData& td
)
{
    //- Call base class function
    CloudType::move(cloud,td,mesh_.time().deltaTValue());
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::checkMassConsistency(scalar& fvM, scalar& ppM)
{
    forAll(mesh_.cells(), celli)
    {
        scalar celliMass = mesh_.V()[celli] * rho_[celli];
        
        fvM += celliMass;
    }

    forAllIters(*this, iter)
    {
        particleType& p = iter();
        
        ppM += p.m();
    }
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::addNewParticle
(
//    const vector& position,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const scalar wt,
    const scalar mass,
    const label patchI,
    const label patchFace,
    const bool iniRls
)
{
    particleType* ptr = new particleType(mesh_,coordinates,cellI,tetFaceI,tetPtI);

    this->addParticle(ptr);
    
    setParticleProperties(*ptr, mass, wt, patchI, patchFace,iniRls);
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    particle.m() = mass;

    particle.wt() = wt;

    particle.sCell() = pManager_.getSuperCellID(particle.cell());
    
    particle.initStatisticalSampling();
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::setEulerianStatistics()
{
    // Do nothing
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::updateEulerianStatistics()
{
    // Do nothing
}


template<class CloudType>
void Foam::ItoPopeCloud<CloudType>::writeFields() const
{
    if (this->size())
    {
        particleType::itoParticleIOType::writeFields(*this);
    }
}




// ************************************************************************* //
