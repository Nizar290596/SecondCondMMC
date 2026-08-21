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
#include "ThermoPopeCloud.H"

#include "CompositionModel.H"
#include "ThermoPhysicalCouplingModel.H"//!!!
#include "interpolationLookUpTable.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::setModels(const mmcVarSet& Xi)
{
    compositionModel_.reset
    (
        CompositionModel<ThermoPopeCloud<CloudType,ReactionThermo> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    thermoPhysicalCouplingModel_.reset
    (
        ThermoPhysicalCouplingModel<ThermoPopeCloud<CloudType,ReactionThermo> >::New
        (
            this->subModelProperties(),
            *this,
            Xi //Since coupling variables are inside submodel
        ).ptr()
    );

}


template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::cloudReset(ThermoPopeCloud<CloudType,ReactionThermo>& c)
{
    CloudType::cloudReset(c);

    compositionModel_ = c.compositionModel_;
    
    thermoPhysicalCouplingModel_ = c.thermoPhysicalCouplingModel_;
}


template<class CloudType, class ReactionThermo>
ReactionThermo*
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::lookupOrConstructThermo()
{
    const fvMesh& mesh = this->mesh();

    const word ReactionThermoName = ReactionThermo::typeName;

    // Find all available ReactionThermo objects in the databank
    auto options = mesh.objectRegistry::lookupClass<ReactionThermo>();

    // Unfortunately this returns a const list of pointers but we need 
    // a non-const access.
    // Solution: Store the name and lookup the pointer with getObjectPtr()

    ReactionThermo* ptr = nullptr;

    // Lookup the name of the thermo model from the cloud submodel dictionary
    const word thermoModelName = 
        this->subModelProperties().template lookupOrDefault<word>
        (
            "thermoModel",
            "thermophysicalProperties"
        );

    // Check if the model can be found in the object registry
    if (mesh.foundObject<ReactionThermo>(thermoModelName))
    {
        Info << "Found thermo model " << thermoModelName 
                 << " in the object registry" << endl;
        ptr = mesh.objectRegistry::getObjectPtr<ReactionThermo>(thermoModelName);
    }
    else
    {
        WarningInFunction 
            << "ThermoPopeCloud requires " << thermoModelName 
            << " thermo model of type " << ReactionThermoName << nl
            << "But model could not be found in objectRegistry. "
            << "Following options are available " << options << nl
            << "You can specify the thermo model name in the cloudProperties "
            << "submodel dictionary with the keyword thermoModel. E.g.," << nl
            << "{"<<nl
            << "    thermoModel   thermophysicalPropertiesCloud;" << nl
            << "    ..." << nl
            << "}" << endl;

        Info << "Attempt to create thermo model "<< thermoModelName << endl;
        
        auto aPtr = ReactionThermo::New(mesh,word::null,thermoModelName);
        ptr = aPtr.release();
        createdThermoModel_ = true;
    }

    return ptr;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, class ReactionThermo>
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::ThermoPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
    :
    CloudType
    (
        cloudName,
        mesh,
        U,
        DEff,
        rho,
        gradRho,
        false,
        true
    ),

    thermoPopeCloud(),

    cloudCopyPtr_(nullptr),

    thermo_
    (
        lookupOrConstructThermo()
    ),

    SLGThermo_(new SLGThermo(mesh,thermo())),

    T_(thermo().T()),

    p_(thermo().p()),

    compositionModel_(nullptr),

    thermoPhysicalCouplingModel_(nullptr)
{
    Info << nl << "Creating thermo Pope Particle Cloud." << nl << endl;

    setModels(Xi);
    
    Info << nl << "Composition and coupling models constructed." << endl;

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Thermo Pope particle cloud data from file." << endl;

            particleType::thermoParticleIOType::readFields(*this, this->composition());
        }
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of thermo Pope particles into the finite volume field." << nl << endl;

                this->initReleaseParticles();
            }
        }
    }
}


template<class CloudType, class ReactionThermo>
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::ThermoPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
    :
    CloudType
    (
        cloudName,
        mesh,
        false,          // Only top level cloud calls initAtCnstr
        readFields
    ),

    thermoPopeCloud(),

    cloudCopyPtr_(nullptr),

    thermo_
    (
        lookupOrConstructThermo()
    ),

    SLGThermo_(new SLGThermo(mesh,thermo())),

    T_(thermo().T()),

    p_(thermo().p()),

    compositionModel_(nullptr),

    thermoPhysicalCouplingModel_(nullptr)
{
    Info << nl << "Creating thermo Pope Particle Cloud." << nl << endl;

    setModels(Xi);
    
    Info << nl << "Composition and coupling models constructed." << endl;

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Thermo Pope particle cloud data from file." << endl;

            particleType::thermoParticleIOType::readFields(*this, this->composition());
        }
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of thermo Pope particles into the finite volume field." << nl << endl;

                this->initReleaseParticles();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, class ReactionThermo>
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::~ThermoPopeCloud()
{
    if (createdThermoModel_)
        delete thermo_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle,mass,wt,patchI,patchFace,iniRls);

    // Set static properties if unset:
    if 
    (
        particleType::indexInXiC_.empty() 
     && particleType::XiCNames_.empty()
     && particleType::componentNames_.empty()
    )
    {
        // ThermoPopeCloud can access the private varaibles as it is a
        // friend class
        particleType::indexInXiC_ = coupling().XiC().cVarInXiC();
        particleType::XiCNames_ = coupling().XiCNames();
        particleType::componentNames_ = composition().componentNames();
    }


    if(!iniRls)
    {
        particle.Y() = composition().YMixture0(patchI,patchFace);

        particle.dpYsource().setSize(particle.Y().size(),0);

        particle.XiC() = coupling().XiC0(patchI,patchFace);

        particle.dpXiCsource().setSize(particle.XiC().size(),0);

        particle.T() = this->T().boundaryField()[patchI][patchFace];

        particle.pc() = this->p().boundaryField()[patchI][patchFace];

    }
    else
    {
        particle.Y() = composition().YMixture0(particle.cell());

        particle.dpYsource().setSize(particle.Y().size(),0);

        particle.XiC() = coupling().XiC0(particle.cell());

        particle.dpXiCsource().setSize(particle.XiC().size(),0);

        particle.T() = this->T()[particle.cell()];

        particle.pc() = this->p()[particle.cell()];

    }

    particle.NumActSp() = 0;

    particle.hA()   = composition().particleMixture(particle.Y()).Ha(particle.pc(),particle.T());

    particle.hEqv() = composition().particleMixture(particle.Y()).Hs(this->p()[particle.cell()],particle.T());

    particle.initStatisticalSampling();
}


template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::setEulerianStatistics()
{

    forAll(this->coupling().XiCNames(),XiCI)
    {
        word XiCName = coupling().XiCNames()[XiCI];

        if (this->eulerianStatsDict().found(XiCName))
        {
            const dimensionSet dim = dimless;// Are all dimless?

            this->eulerianStats().newProperty(XiCName,dim);
        }
    }

    if (this->eulerianStatsDict().found("T"))
    {
        const dimensionSet dim = dimTemperature;

        this->eulerianStats().newProperty("T",dim);
    }

    forAll(this->composition().componentNames(),specieI)
    {
        word specieName = composition().componentNames()[specieI];

        if (this->eulerianStatsDict().found(specieName))
        {
            const dimensionSet dim = dimless;

            this->eulerianStats().newProperty("Y"+specieName,dim);

            if (composition().printMoleFractionsEnabled())
                this->eulerianStats().newProperty("X"+specieName,dim);
        }
    }
}


template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();

    forAllIters(*this, iter)
    {
        if (this->secondCondMixingEnabled() && iter().secondCondFlag() != 1)
            continue;
	this->eulerianStats().findCell(iter().position());

        if (this->eulerianStatsDict().found("T"))
            this->eulerianStats().calculate("T",iter().wt(),iter().T());

        forAll(composition().componentNames(),specieI)
        {
            const word& specieName = composition().componentNames()[specieI];

            if (this->eulerianStatsDict().found(specieName))
            {
                this->eulerianStats().calculate
                (
                    "Y" +specieName,iter().wt(),
                    iter().Y()[specieI]
                );

                if (!composition().printMoleFractionsEnabled())
                    continue;

                scalar moleFraction = composition().X(iter().Y())[specieI];

                this->eulerianStats().calculate
                (
                    "X" + specieName,iter().wt(),
                    moleFraction
                );
             }
        }

        //- statistics of coupling variables
        forAll(this->coupling().XiCNames(),XiCI)
        {
            const word& couplingName = coupling().XiCNames()[XiCI];

            if (this->eulerianStatsDict().found(couplingName))
                this->eulerianStats().calculate(couplingName,iter().wt(),iter().XiC(couplingName));
        }
    }
}


template<class CloudType, class ReactionThermo>
void Foam::ThermoPopeCloud<CloudType,ReactionThermo>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::thermoParticleIOType::writeFields(*this, this->composition());
    }
}

