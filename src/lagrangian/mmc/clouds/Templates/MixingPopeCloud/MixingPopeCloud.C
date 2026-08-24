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

#include "MixingPopeCloud.H"
#include "CloudMixingModel.H"
#include "interpolationCellPoint.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setModels(const mmcVarSet& Xi)
{
    mixingModel_.reset
    (
        CloudMixingModel<MixingPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this,
            Xi //Since reference variables are inside submodel!!!
        ).ptr()
    );

    // Conditionally construct the second-conditioning mixing model.
    // The secondConditioning sub-dictionary must be present in cloudProperties
    // and have   enabled true;   for the model to be activated.
    if (this->cloudProperties().found("secondConditioning"))
    {
        const dictionary& scDict =
            this->cloudProperties().subDict("secondConditioning");
 
        if (scDict.lookupOrDefault("enabled", false))
        {
            const word scModelType
            (
                this->subModelProperties().lookup("secondCondMixingModel")
            );
 
            auto cstrIter =
                CloudMixingModel<MixingPopeCloud<CloudType>>::
                    dictionaryConstructorTablePtr_->find(scModelType);
 
            if (cstrIter == CloudMixingModel<MixingPopeCloud<CloudType>>::dictionaryConstructorTablePtr_->end())
            {
                FatalErrorInFunction
                    << "Unknown secondCondMixingModel type "
                    << scModelType << nl
                    << "Valid types are:" << nl
                    << CloudMixingModel<MixingPopeCloud<CloudType>>::
                           dictionaryConstructorTablePtr_->sortedToc()
                    << exit(FatalError);
            }
 
            secondCondMixingModel_.reset
            (
                cstrIter()(this->subModelProperties(), *this, Xi)
            );
 
            Info << "Second-conditioning mixing model: "
                 << scModelType << endl;

            // The W(phi) = A*(1-phi)*exp[Z*(phi-1)] progress-variable source
            // has been removed: it was non-negative for every phi <= 1, so it
            // drove all particles to phi = 1 regardless of local conditions.
            // phi is now relaxed towards the Eulerian progress variable c with
            // timescale tauC instead. Fail loudly on cases that still set the
            // old keys rather than silently changing their behaviour.
            if
            (
                scDict.found("A_phi")
             || scDict.found("Z_phi")
            )
            {
                FatalErrorInFunction
                    << "secondConditioning contains A_phi/Z_phi, but the "
                    << "W(phi) progress-variable source has been removed." << nl
                    << "phi is now relaxed towards the Eulerian progress "
                    << "variable c over a timescale tauC." << nl
                    << "Remove A_phi and Z_phi and set tauC instead."
                    << exit(FatalError);
            }

            if (!scDict.found("tauC"))
            {
                WarningInFunction
                    << "secondConditioning does not set tauC; defaulting to 0, "
                    << "i.e. phi is set equal to the Eulerian progress "
                    << "variable c at the particle position every time step."
                    << endl;
            }
        }
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::cloudReset(MixingPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
    
    mixingModel_.reset(c.mixingModel_.ptr());

    secondCondMixingModel_.reset(c.secondCondMixingModel_.ptr());
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixingPopeCloud<CloudType>::MixingPopeCloud
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
        Xi,
        false,  // Only the top level cloud will call initAtCnstr
        readFields
    ),
  
    mixingPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    mixingModel_(nullptr),

	secondCondMixingModel_(nullptr),

    secondCondR_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("R", 0.0)
    ),
    secondCondBeta_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("beta", 1.0)
    ),
    secondCondTauOU_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("tauOU", 1.0)
    ),
    secondCondTauC_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("tauC", 0.0)
    ),
    secondCondTu_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Tu", 300.0)
    ),
    secondCondTb_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Tb", 2000.0)
    ),
    phiDeviation_(0.0)

{
    Info << "Creating mixing Pope Particle Cloud." << nl << endl;

    setModels(Xi); // passing mmcVarSet
    
    Info << nl << "Mixing model constructed." << endl;
    
    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Mixing Pope particle cloud data from file." << endl;
            
            particleType::mixingParticleIOType::readFields(*this, this->mixing());
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial realease of Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
            }
        }
    }
}


template<class CloudType>
Foam::MixingPopeCloud<CloudType>::MixingPopeCloud
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
        Xi,
        false,      // Only the top level cloud will call initAtCnstr
        readFields
    ),
  
    mixingPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    mixingModel_(nullptr),
    
    secondCondMixingModel_(nullptr),

    secondCondR_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("R", 0.0)
    ),
    secondCondBeta_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("beta", 1.0)
    ),
    secondCondTauOU_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("tauOU", 1.0)
    ),
    secondCondTauC_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("tauC", 0.0)
    ),
    secondCondTu_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Tu", 300.0)
    ),
    secondCondTb_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Tb", 2000.0)
    ),
    phiDeviation_(0.0)
    
{
    Info << "Creating mixing Pope Particle Cloud." << nl << endl;

    setModels(Xi); // passing mmcVarSet
    
    Info << nl << "Mixing model constructed." << endl;
    
    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Mixing Pope particle cloud data from file." << endl;
            
            particleType::mixingParticleIOType::readFields(*this, this->mixing());
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial realease of Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MixingPopeCloud<CloudType>::~MixingPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::updatePhi(const scalar deltaT)
{
    // Relax the particle progress variable towards the Eulerian progress
    // variable c interpolated at the particle position:
    //
    //     phi <- phi + (1 - exp(-dt/tauC)) * (c_p - phi)
    //
    // The update is the exact solution of dphi/dt = (c_p - phi)/tauC for a
    // constant c_p, so it carries no truncation error and is independent of
    // dt -- the same construction used by OUStateUpdate() for omega.
    //
    // Applied to ALL particles (not only flagged ones) so phi is consistent
    // across the whole cloud before the second conditioning runs, and applied
    // AFTER mixing().Smix() so it acts on the mixed phi.
    //
    // This replaces the former W(phi) = A*(1-phi)*exp[Z*(phi-1)] source. That
    // term was non-negative for every phi <= 1 and had no dependence on
    // temperature, mixture or proximity to a flame, so it drove every particle
    // in the domain monotonically to phi = 1 -- making phi a function of
    // residence time rather than of reaction, and leaving
    // phi_deg = phi*exp(beta*omega) as pure OU noise once it saturated.
    const scalar tauC = secondCondTauC_;

    if (tauC < 0)
        return;

    // tauC == 0 is the limiting case of instantaneous relaxation, i.e. phi is
    // a pure interpolation of c (the contract of referenceType 'interpolated'
    // in mmcPremixedFoam). Guard the exponential against division by zero.
    const scalar relax =
        (tauC > SMALL) ? (1.0 - Foam::exp(-deltaT/tauC)) : 1.0;

    // The Eulerian progress variable is registered by the solver's
    // createFields.H. It is deliberately not part of the mmcVarSet, so it is
    // looked up from the object registry here rather than through Xi.
    if (!this->mesh().objectRegistry::foundObject<volScalarField>("c"))
    {
        FatalErrorInFunction
            << "Second conditioning is enabled but the Eulerian progress "
            << "variable field 'c' is not registered." << nl
            << "The solver must create it (see SPFoam createFields.H) and "
            << "update it each time step before the particles are evolved."
            << exit(FatalError);
    }

    const volScalarField& cField =
        this->mesh().objectRegistry::lookupObject<volScalarField>("c");

    interpolationCellPoint<scalar> c_intp_(cField);

    // Diagnostics: mean |phi - c_p| over the cloud, so drift of phi away from
    // the resolved progress variable is directly observable in the log.
    scalar sumWt = 0.0;
    scalar sumAbsDev = 0.0;

    forAllIters(*this, iter)
    {
        const scalar cP = c_intp_.interpolate
        (
            iter().position(), iter().cell(), iter().face()
        );

        scalar& phi = iter().phi();

        phi += relax*(cP - phi);
        phi  = max(0.0, min(1.0, phi));

        const scalar w = iter().wt();
        sumWt     += w;
        sumAbsDev += w*mag(phi - cP);

        // For non-subset particles (omegaOU == 0 always), phiModified = phi.
        // For subset particles, updateOUProcess() recomputes phiModified as
        // phi*exp(beta*omegaOU) immediately after this call.
        if (iter().secondCondFlag() != 1)
            iter().phiModified() = phi;
    }

    reduce(sumWt, sumOp<scalar>());
    reduce(sumAbsDev, sumOp<scalar>());

    phiDeviation_ = (sumWt > SMALL) ? sumAbsDev/sumWt : 0.0;

    Info<< "    [secondCond] mean |phi - c| = " << phiDeviation_
        << " (tauC = " << tauC << " s, relax = " << relax << ")" << endl;
}
 
 
template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::updateOUProcess(const scalar deltaT)
{
    // For each particle flagged for second conditioning:
    //   1. Advance ω_OU using the exact discrete OU update.
    //   2. Recompute φ° = φ·exp(β·ω_OU) so buildParticleList() in
    //      secondCondMixing().Smix() sees the current modified variable.
    const scalar beta  = secondCondBeta_;
    const scalar tauOU = secondCondTauOU_;
 
    if (tauOU <= SMALL)
        return;
 
    forAllIters(*this, iter)
    {
        if (iter().secondCondFlag() == 1)
        {
            const scalar xi = this->rndGen_.Normal(0, 1);
            iter().omegaOU() = OUStateUpdate
            (
                iter().omegaOU(), deltaT, tauOU, xi
            );
            iter().phiModified() =
                iter().phi() * Foam::exp(beta * iter().omegaOU());
        }
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle, mass, wt, patchI, patchFace,iniRls);    
    
    /*if 
    (
        particleType::indexInXiR_.empty()
     && particleType::XiRNames_.empty()
    )
    {
        particleType::indexInXiR_ = mixing().XiR().rVarInXiR();
        particleType::XiRNames_   = mixing().XiRNames();
    }*/


    particle.dx() = 0;


    // Assign second-conditioning subset flag based on fraction R
    particle.secondCondFlag() =
        (secondCondR_ > 0 && this->rndGen_.Random() < secondCondR_) ? 1 : 0;
 
    // Initialize progress variable from particle temperature
    // phi = (T - Tu) / (Tb - Tu), clamped to [0, 1].
    // This is only the initial condition: from here on phi is relaxed towards
    // the Eulerian progress variable c by updatePhi() and mixed pairwise by
    // the first-conditioning model.
    {
        const scalar dT = secondCondTb_ - secondCondTu_;
        if (dT > SMALL)
        {
            particle.phi() =
                max(0.0, min(1.0, (particle.T() - secondCondTu_) / dT));
        }
        else
        {
            particle.phi() = 0.0;
        }
    }
 
//    particle.phiModified() = particle.phi();

    // Initialise omegaOU from the OU stationary distribution N(0,1) for
    // flagged (subset) particles only; non-subset particles never get
    // OU-updated downstream so their omegaOU stays at 0.
    if (particle.secondCondFlag() == 1)
    {
        particle.omegaOU() = this->rndGen_.Normal(0, 1);
    }
 
    particle.phiModified() =
        particle.phi() * Foam::exp(secondCondBeta_ * particle.omegaOU());


    label numXiR = mixing().numXiR();

    //- Extension for a set of variables    
    particle.dXiR().setSize(numXiR,0.0);
    
    if (!iniRls)
    {
        particle.XiR() = mixing().XiR0(patchI,patchFace,particle);
    }
    else
    {
        particle.XiR() = mixing().XiR0(particle.cell(),particle);
    } 

    particle.initStatisticalSampling();
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setEulerianStatistics()
{
    
    if (this->eulerianStatsDict().found("dx"))
    {
        const dimensionSet dim = dimLength;
        
        this->eulerianStats().newProperty("dx",dim);
    }
    
    //- statistics of reference variables (mixing distances)
    forAll(this->mixing().XiRNames(),XiRI)
    {
        word refVarName = "d" + mixing().XiRNames()[XiRI];
        
        if (this->eulerianStatsDict().found(refVarName))
        {
            const dimensionSet dim = dimless;// they will be considered dimless 
            
            this->eulerianStats().newProperty(refVarName,dim);   
        }
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();

    forAllIters(*this, iter)
    {
        this->eulerianStats().findCell(iter().position());

        //- statistics of distance for reference variables 
        forAll(mixing().XiRNames(),XiRI)
        {
            const word& refVarName = mixing().XiRNames()[XiRI];
            const word& dRefVarName = "d" + refVarName;
            
            if (this->eulerianStatsDict().found(dRefVarName))
                this->eulerianStats().calculate(dRefVarName,iter().wt(),iter().dXiR(refVarName));
        }
            
        if (this->eulerianStatsDict().found("dx"))
            this->eulerianStats().calculate("dx",iter().wt(),iter().dx());
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::mixingParticleIOType::writeFields(*this, this->mixing());
    }
}

// ************************************************************************* //

