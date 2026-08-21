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

#include "SprayKinematicCloud.H"
#include "interpolation.H"
#include "subCycleTime.H"

#include "InjectionModelList.H"
#include "DispersionModel.H"
#include "PatchInteractionModel.H"
#include "StochasticCollisionModel.H"
#include "SurfaceFilmModel.H"
#include "profiling.H"

#include "PackingModel.H"
#include "ParticleStressModel.H"
#include "DampingModel.H"
#include "IsotropyModel.H"
#include "TimeScaleModel.H"

#include "surfaceFields.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::setModels()
{
    dispersionModel_.reset
    (
        DispersionModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    patchInteractionModel_.reset
    (
        PatchInteractionModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    stochasticCollisionModel_.reset
    (
        StochasticCollisionModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    surfaceFilmModel_.reset
    (
        SurfaceFilmModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    packingModel_.reset
    (
        PackingModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    dampingModel_.reset
    (
        DampingModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    isotropyModel_.reset
    (
        IsotropyModel<SprayKinematicCloud<CloudType>>::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
 
    UIntegrator_.reset
    (
        integrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
        ).ptr()
    );
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::parcelControl()
{
    Info << "kill parcels for extremely small diamter" <<endl;
           
    forAllIters(*this, iter)
    {
        parcelType& p = iter();      
                   
        if ( p.d() < dMin_)
        {
            this->deleteParticle(iter());   //- delete the parcels with extremely small diameter     
        }
    }
} 

template<class CloudType>
template<class TrackCloudType>
void Foam::SprayKinematicCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    if (solution_.steadyState())
    {
        cloud.storeState();

        cloud.preEvolve(td);

        evolveCloud(cloud,td);

        if (solution_.coupled())
        {
            cloud.relaxSources(cloud.cloudCopy());
        }
    }
    else
    {
        cloud.preEvolve(td);

        evolveCloud(cloud,td);

        if (solution_.coupled())
        {
            cloud.scaleSources();
        }
    }

    cloud.info();

    cloud.postEvolve(td);

    if (solution_.steadyState())
    {
        cloud.restoreState();
    }

    //- parcel management   
    this->parcelControl(); 
   
    //-Eulerian statistics update    
    this->updateEulerianStatistics();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::buildCellOccupancy()
{
    if (!cellOccupancyPtr_.good())
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<parcelType*>>(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size

        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    List<DynamicList<parcelType*>>& cellOccupancy = cellOccupancyPtr_();

    forAll(cellOccupancy, cO)
    {
        cellOccupancy[cO].clear();
    }

    forAllIters(*this, iter)
    {
        cellOccupancy[iter().cell()].append(&iter());
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::updateCellOccupancy()
{
    // Only build the cellOccupancy if the pointer is set, i.e. it has
    // been requested before.
    if (cellOccupancyPtr_.valid())
    {
        buildCellOccupancy();
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::SprayKinematicCloud<CloudType>::evolveCloud
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    if (solution_.coupled())
    {
        cloud.resetSourceTerms();
    }

    if (solution_.transient())
    {
        label preInjectionSize = this->size();

//        this->surfaceFilm().inject(td);

        // Update the cellOccupancy if the size of the cloud has changed
        // during the injection.
        if (preInjectionSize != this->size())
        {
            
            updateCellOccupancy();
            
            preInjectionSize = this->size();
        }

        injectors_.inject(cloud,td);


        // Assume that motion will update the cellOccupancy as necessary
        // before it is required.
        cloud.motion(cloud,td);
    }
    else
    {
//        this->surfaceFilm().injectSteadyState(td);

        injectors_.injectSteadyState(cloud, td, solution_.trackTime());

        td.part() = parcelType::trackingData::tpLinearTrack;
        CloudType::move(cloud, td,  solution_.trackTime());
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    Info<< endl;

    if (debug)
    {
        this->writePositions();
    }

    forces_.cacheFields(false);

    functions_.postEvolve(td);

    solution_.nextIter();

    if (this->db().time().writeTime())
    {
        outputProperties_.regIOobject::write();
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::cloudReset(SprayKinematicCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    rndGen_ = c.rndGen_;

    forces_.transfer(c.forces_);

    functions_.transfer(c.functions_);

    injectors_.transfer(c.injectors_);

    patchInteractionModel_.reset(c.patchInteractionModel_.ptr());

    UIntegrator_.reset(c.UIntegrator_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayKinematicCloud<CloudType>::SprayKinematicCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& muEff,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(rho.mesh(), cloudName, false),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fuelCloud1Properties_
    (
        IOobject
        (
            cloudName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    constProps_(particleProperties_, solution_.active()),
    subModelProperties_
    (
        particleProperties_.subOrEmptyDict("subModels", keyType::REGEX, solution_.active())
    ),
    rndGen_(Pstream::myProcNo()),
    cellOccupancyPtr_(),
    cellLengthScale_(mag(cbrt(mesh_.V()))),
    rho_(rho),
    U_(U),
    mu_(muEff),
    g_(g),
    pAmbient_(0.0),
    forces_
    (
        *this,
        mesh_,
        subModelProperties_.subOrEmptyDict
        (
            "particleForces",
            keyType::REGEX,
            solution_.active()
        ),
        solution_.active()
    ),
    functions_
    (
        *this,
        particleProperties_.subOrEmptyDict("cloudFunctions"),
        solution_.active()
    ),
    injectors_
    (
        subModelProperties_.subOrEmptyDict("injectionModels"),
        *this
    ),
    dispersionModel_(nullptr),
    patchInteractionModel_(nullptr),
    stochasticCollisionModel_(nullptr),
    surfaceFilmModel_(nullptr),
    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr),
    UIntegrator_(nullptr),
    UTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + "UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, Zero)
        )
    ),
    UCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimMass, 0.0)
        )
    ),
    parcelManagementDict_
    (
        fuelCloud1Properties_.subDict("parcelManagement")
    ),
    dMin_(readScalar(parcelManagementDict_.lookup("dMin"))),
    eulerianStatsDict_
    (
        particleProperties_.subDict("eulerianStatistics")
    ),
    eulerianStats_(mesh_,eulerianStatsDict_)
{
    if (solution_.active())
    {
        setModels();

        setEulerianStatistics();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (solution_.resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::SprayKinematicCloud<CloudType>::SprayKinematicCloud
(
    SprayKinematicCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c.mesh_, name, c),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(c.mesh_),
    particleProperties_(c.particleProperties_),
    fuelCloud1Properties_(c.fuelCloud1Properties_),
    outputProperties_(c.outputProperties_),
    solution_(c.solution_),
    constProps_(c.constProps_),
    subModelProperties_(c.subModelProperties_),
    rndGen_(c.rndGen_, true),
    cellOccupancyPtr_(nullptr),
    cellLengthScale_(c.cellLengthScale_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu()),
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    forces_(c.forces_),
    functions_(c.functions_),
    injectors_(c.injectors_),
    dispersionModel_(c.dispersionModel_->clone()),
    patchInteractionModel_(c.patchInteractionModel_->clone()),
    stochasticCollisionModel_(c.stochasticCollisionModel_->clone()),
    surfaceFilmModel_(c.surfaceFilmModel_->clone()),
    packingModel_(c.packingModel_->clone()),
    dampingModel_(c.dampingModel_->clone()),
    isotropyModel_(c.isotropyModel_->clone()),
    UIntegrator_(c.UIntegrator_->clone()),
    UTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + "UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.UTrans_()
        )
    ),
    UCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                name + "UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.UCoeff_()
        )
    ),
    parcelManagementDict_
    (
        fuelCloud1Properties_.subDict("parcelManagement")
    ),
    dMin_(c.dMin_),
    eulerianStatsDict_
    (
        particleProperties_.subDict("eulerianStatistics")
    ),
    eulerianStats_(mesh_,eulerianStatsDict_)
{}


template<class CloudType>
Foam::SprayKinematicCloud<CloudType>::SprayKinematicCloud
(
    const fvMesh& mesh,
    const word& name,
    const SprayKinematicCloud<CloudType>& c
)
:
    CloudType(mesh, name, IDLList<parcelType>()),
    kinematicCloud(),
    cloudCopyPtr_(nullptr),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    fuelCloud1Properties_
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    outputProperties_
    (
        IOobject
        (
            name + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/name,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    solution_(mesh),
    constProps_(),
    subModelProperties_(dictionary::null),
    rndGen_(),
    cellOccupancyPtr_(nullptr),
    cellLengthScale_(c.cellLengthScale_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu()),
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    forces_(*this, mesh),
    functions_(*this),
    injectors_(*this),
    dispersionModel_(nullptr),
    patchInteractionModel_(nullptr),
    stochasticCollisionModel_(nullptr),
    surfaceFilmModel_(nullptr),
    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr),
    UIntegrator_(nullptr),
    UTrans_(nullptr),
    UCoeff_(nullptr),
    parcelManagementDict_
    (
        fuelCloud1Properties_.subDict("parcelManagement")
    ), 
    dMin_(c.dMin_),
    eulerianStatsDict_
    (
        particleProperties_.subDict("eulerianStatistics")
    ),
    eulerianStats_(mesh_, eulerianStatsDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayKinematicCloud<CloudType>::~SprayKinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::SprayKinematicCloud<CloudType>::hasWallImpactDistance() const
{
    return true;
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    parcel.rho() = constProps_.rho0();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    const scalar carrierDt = mesh_.time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;

    if (parcel.typeId() == -1)
    {
        parcel.typeId() = constProps_.parcelTypeId();
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<SprayKinematicCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::resetSourceTerms()
{
    UTrans().field() = Zero;
    UCoeff().field() = 0.0;
}


template<class CloudType>
template<class Type>
void Foam::SprayKinematicCloud<CloudType>::relax
(
    DimensionedField<Type, volMesh>& field,
    const DimensionedField<Type, volMesh>& field0,
    const word& name
) const
{
    const scalar coeff = solution_.relaxCoeff(name);
    field = field0 + coeff*(field - field0);
}


template<class CloudType>
template<class Type>
void Foam::SprayKinematicCloud<CloudType>::scale
(
    DimensionedField<Type, volMesh>& field,
    const word& name
) const
{
    const scalar coeff = solution_.relaxCoeff(name);
    field *= coeff;
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::relaxSources
(
    const SprayKinematicCloud<CloudType>& cloudOldTime
)
{
    this->relax(UTrans_(), cloudOldTime.UTrans(), "U");
    this->relax(UCoeff_(), cloudOldTime.UCoeff(), "U");
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::scaleSources()
{
    this->scale(UTrans_(), "U");
    this->scale(UCoeff_(), "U");
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    // force calculaion of mesh dimensions - needed for parallel runs
    // with topology change due to lazy evaluation of valid mesh dimensions
    label nGeometricD = mesh_.nGeometricD();

    Info<< "\nSolving " << nGeometricD << "-D cloud " << this->name() << endl;

    forces_.cacheFields(true);
    updateCellOccupancy();

    pAmbient_ = constProps_.dict().template
        lookupOrDefault<scalar>("pAmbient", pAmbient_);

    functions_.preEvolve(td);
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::evolve()
{
    if (solution_.canEvolve())
    {
        typename parcelType::trackingData td(*this);
        solve(*this,td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::SprayKinematicCloud<CloudType>::motion
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)

{
    td.part() = parcelType::trackingData::tpLinearTrack;
    CloudType::move(cloud, td, solution_.trackTime());

    updateCellOccupancy();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::patchData
(
    const parcelType& p,
    const polyPatch& pp,
    vector& nw,
    vector& Up
) const
{
    p.patchData(nw, Up);

    // If this is a wall patch, then there may be a non-zero tangential velocity
    // component; the lid velocity in a lid-driven cavity case, for example. We
    // want the particle to interact with this velocity, so we look it up in the
    // velocity field and use it to set the wall-tangential component.
    if (isA<wallPolyPatch>(pp))
    {
        const label patchi = pp.index();
        const label patchFacei = pp.whichFace(p.face());

        // We only want to use the boundary condition value  onlyif it is set
        // by the boundary condition. If the boundary values are extrapolated
        // (e.g., slip conditions) then they represent the motion of the fluid
        // just inside the domain rather than that of the wall itself.
        if (U_.boundaryField()[patchi].fixesValue())
        {
            const vector Uw1 = U_.boundaryField()[patchi][patchFacei];
            const vector& Uw0 =
                U_.oldTime().boundaryField()[patchi][patchFacei];

            const scalar f = p.currentTimeFraction();

            const vector Uw = Uw0 + f*(Uw1 - Uw0);

            const tensor nnw = nw*nw;

            Up = (nnw & Up) + Uw - (nnw & Uw);
        }
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::updateMesh()
{
    updateCellOccupancy();
    injectors_.updateMesh();
    cellLengthScale_ = mag(cbrt(mesh_.V()));
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    updateMesh();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::info()
{
    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    Info<< "Cloud: " << this->name() << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Linear momentum                 = "
        << linearMomentum << nl
        << "   |Linear momentum|                = "
        << mag(linearMomentum) << nl
        << "    Linear kinetic energy           = "
        << linearKineticEnergy << nl;

    injectors_.info();
    this->patchInteraction().info();
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::setEulerianStatistics()
{
    for (label i = 0; i < 6; i++)
    {
        word Ud_xyz;
               
        if (i == 0)
            Ud_xyz = "u";          

        else if (i == 1)
            Ud_xyz = "v";           
            
        else if (i == 2)
            Ud_xyz = "w";

        else if (i == 3)
            Ud_xyz = "uAll";    
        
        else if (i == 4)
            Ud_xyz = "vAll";

        else if (i == 5)
            Ud_xyz = "wAll";

        if (this->eulerianStatsDict().found(Ud_xyz))
        {
            const dimensionSet dim = dimVelocity;
        
            typedef VectorSpace<Vector<scalar>, scalar, 2> pair;
        

            List<pair> d1d2(this->eulerianStatsDict().lookup(Ud_xyz));
       
            forAll(d1d2, nd)
            {
                const word rangeName(name(d1d2[nd][0]*1e6) + "to" + name(d1d2[nd][1]*1e6) );

                this->eulerianStats().newProperty(Ud_xyz+rangeName,dim);
            }         
        }       
                
        
    }
}


template<class CloudType>
void Foam::SprayKinematicCloud<CloudType>::updateEulerianStatistics()
{ 
    for (label i = 0; i < 6; i++)
    {
        word Ud_xyz;
                
        if (i == 0)
            Ud_xyz = "u";           

        else if (i == 1)
            Ud_xyz = "v";            
            
        else if (i == 2)
            Ud_xyz = "w";

        else if (i == 3)
            Ud_xyz = "uAll";    
        
        else if (i == 4)
            Ud_xyz = "vAll";

        else if (i == 5)
            Ud_xyz = "wAll";

        if (this->eulerianStatsDict().found(Ud_xyz))
        {
            typedef VectorSpace<Vector<scalar>, scalar, 2> pair;
            
            List<pair> d1d2(this->eulerianStatsDict().lookup(Ud_xyz));
        
            forAllIters(*this, iter)
            {
                forAll(d1d2, nd)
                {
                    if (iter().d() >= d1d2[nd][0] && iter().d() < d1d2[nd][1])
                    {
                        const word rangeName(name(d1d2[nd][0]*1e6) + "to" + name(d1d2[nd][1]*1e6) );
                        
                        scalar vel;
                        if (i == 0)
                            vel = iter().U().x();
                        else if (i == 1)
                            vel = iter().U().y();
                        else if (i == 2)
                            vel = iter().U().z();
                        else if (i == 3) // Ux-All Sizes
                            vel = iter().U().x();
                        else if (i == 4) // Uy-All Sizes
                            vel = iter().U().y();
                        else if (i == 5) // Uz-All Sizes
                            vel = iter().U().z();
            
                       this->eulerianStats().findCell(iter().position());

                        scalar wt = 1.0;
                        
                        this->eulerianStats().calculate(Ud_xyz+rangeName,wt,vel);
                        
                        break;
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
