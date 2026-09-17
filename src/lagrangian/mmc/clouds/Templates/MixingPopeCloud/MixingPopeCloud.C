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


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::setModels(const mmcVarSet& Xi)
{
    const dictionary sc = this->cloudProperties().subOrEmptyDict("secondConditioning");
    if (sc.lookupOrDefault("enabled", false))
    {
        if (sc.found("referenceDiffusionRate"))
        {
            const scalar rate = readScalar(sc.lookup("referenceDiffusionRate"));
            if (sc.found("beta") || !(rate >= 0 && std::isfinite(rate)) || secondCondTauOU_ <= 0)
                FatalErrorInFunction << "Use beta OR referenceDiffusionRate >= 0, with positive tauOU" << exit(FatalError);
            secondCondBeta_ = sqrt(rate*secondCondTauOU_);
        }
        secondCondMaxSourceStep_ = sc.lookupOrDefault<scalar>("maxPhiSourceStep", 0.05);
        secondCondAgeRate_ = sc.lookupOrDefault<scalar>("burnedAgeRate", 0);
        secondCondAgeThreshold_ = sc.lookupOrDefault<scalar>("burnedAgeThreshold", 0.99);
        if (!(secondCondR_ > 0 && secondCondR_ <= 1 && secondCondBeta_ >= 0
            && secondCondTauOU_ >= 0 && (secondCondBeta_ == 0 || secondCondTauOU_ > 0)
            && secondCondTb_ > secondCondTu_ && secondCondAPhi_ >= 0 && secondCondZPhi_ >= 0
            && secondCondMaxSourceStep_ > 0 && secondCondMaxSourceStep_ <= 1
            && secondCondAgeRate_ >= 0 && secondCondAgeThreshold_ >= 0 && secondCondAgeThreshold_ < 1)
            || !std::isfinite(secondCondR_ + secondCondBeta_ + secondCondTauOU_
                + secondCondTb_ + secondCondTu_ + secondCondAPhi_ + secondCondZPhi_
                + secondCondAgeRate_ + secondCondMaxSourceStep_))
            FatalErrorInFunction << "Invalid secondConditioning parameters: require 0<R<=1, beta>=0, "
                << "positive tauOU when beta>0, Tb>Tu, A_phi/Z_phi>=0, "
                << "0<maxPhiSourceStep<=1 and nonnegative aging rate." << exit(FatalError);
        const int seeds[] = {int(sc.lookupOrDefault<label>("randomSeed", 5489)),
                             int(Pstream::myProcNo()), 139691};
        secondCondRandom_.RandomInitByArray(seeds, 3);
        Info<< "Second conditioning: R=" << secondCondR_ << ", A_phi=" << secondCondAPhi_
            << ", Z_phi=" << secondCondZPhi_ << ", beta=" << secondCondBeta_
            << ", tauOU=" << secondCondTauOU_ << ", beta^2/tauOU="
            << (secondCondTauOU_ > 0 ? sqr(secondCondBeta_)/secondCondTauOU_ : 0)
            << ", maxPhiSourceStep=" << secondCondMaxSourceStep_ << nl;
        if (secondCondAgeRate_ > 0)
            WarningInFunction << "burnedAgeRate enables an additional age closure, not an equation "
                << "specified by the paper; validate it against slow-species kinetics." << nl;
    }
    mixingModel_.reset
    (
        CloudMixingModel<MixingPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this,
            Xi //Since reference variables are inside submodel!!!
        ).ptr()
    );

    if (sc.lookupOrDefault("enabled", false)
        && word(this->subModelProperties().lookup("mixingModel")) != "MMCcurl")
        FatalErrorInFunction << "Second conditioning currently requires MMCcurl at the dense level" << exit(FatalError);

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
 
            if (scModelType != "secondCondMMCcurl")
                FatalErrorInFunction << "Supported sparse second-conditioning model is secondCondMMCcurl" << exit(FatalError);

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
        }
    }
}


template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::cloudReset(MixingPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
    
    mixingModel_.reset(c.mixingModel_.ptr());

    secondCondMixingModel_.reset(c.secondCondMixingModel_.ptr());
    secondCondRandom_ = c.secondCondRandom_;
    
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
    secondCondAPhi_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("A_phi", 0.0)
    ),
    secondCondZPhi_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Z_phi", 0.0)
    )

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
    secondCondAPhi_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("A_phi", 0.0)
    ),
    secondCondZPhi_
    (
        this->cloudProperties_.subOrEmptyDict("secondConditioning")
            .template lookupOrDefault<scalar>("Z_phi", 0.0)
    )
    
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
void Foam::MixingPopeCloud<CloudType>::updatePhiReaction(const scalar deltaT)
{
    forAllIters(*this, iter)
    {
        try
        {
            iter().phi() = secondConditioningNumerics::progressReaction
            (
                iter().phi(), deltaT, secondCondAPhi_, secondCondZPhi_, secondCondMaxSourceStep_
            );
        }
        catch (const std::exception& error)
        {
            FatalErrorInFunction << error.what() << exit(FatalError);
        }
    }
}
 
 
template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::updateOUProcess(const scalar deltaT)
{
    forAllIters(*this, iter)
    {
        if (iter().secondCondFlag() == 1)
        {
            if (secondCondBeta_ > 0)
                iter().omegaOU() = OUStateUpdate(iter().omegaOU(), deltaT,
                    secondCondTauOU_, secondCondRandom_.Normal(0, 1));
            else
                iter().omegaOU() = 0;
        }
        // Optional accumulated residence time above a smooth burned-state gate.
        // This is a declared extension, not a calibrated closure from the paper.
        if (secondCondAgeRate_ > 0)
            iter().burnedAge() += deltaT*max(scalar(0),
                (iter().phi()-secondCondAgeThreshold_)/(1-secondCondAgeThreshold_));
        try
        {
            iter().phiModified() = secondConditioningNumerics::modifiedProgress
            (
                iter().phi(), secondCondBeta_, iter().omegaOU(),
                secondCondAgeRate_*iter().burnedAge()
            );
        }
        catch (const std::exception& error)
        {
            FatalErrorInFunction << error.what() << exit(FatalError);
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
    initializeSecondConditioningState(particle);

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
    
    for (const word key : {word("phi"), word("phiModified"), word("omegaOU"), word("secondCondFlag"), word("burnedAge")})
        if (this->eulerianStatsDict().found(key))
            this->eulerianStats().newProperty(key, key == "burnedAge" ? dimTime : dimless);

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
        if (this->eulerianStatsDict().found("phi"))
            this->eulerianStats().calculate("phi",iter().wt(),iter().phi());
        if (this->eulerianStatsDict().found("phiModified"))
            this->eulerianStats().calculate("phiModified",iter().wt(),iter().phiModified());
        if (this->eulerianStatsDict().found("burnedAge"))
            this->eulerianStats().calculate("burnedAge",iter().wt(),iter().burnedAge());
        if (this->eulerianStatsDict().found("secondCondFlag"))
            this->eulerianStats().calculate("secondCondFlag",iter().wt(),scalar(iter().secondCondFlag()));
        if (iter().secondCondFlag() == 1 && this->eulerianStatsDict().found("omegaOU"))
            this->eulerianStats().calculate("omegaOU",iter().wt(),iter().omegaOU());
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



template<class CloudType>
void Foam::MixingPopeCloud<CloudType>::initializeSecondConditioningState(particleType& particle)
{
    particle.secondCondFlag() = 0;
    particle.omegaOU() = 0;
    particle.phi() = 0;
    particle.phiModified() = 0;
    particle.burnedAge() = 0;
    if (!secondCondMixingEnabled()) return;
    particle.secondCondFlag() = secondCondRandom_.Random() < secondCondR_ ? 1 : 0;
    particle.phi() = max(scalar(0), min(scalar(1),
        (particle.T()-secondCondTu_)/(secondCondTb_-secondCondTu_)));
    if (particle.secondCondFlag() == 1 && secondCondBeta_ > 0)
        particle.omegaOU() = secondCondRandom_.Normal(0, 1);
    particle.phiModified() = secondConditioningNumerics::modifiedProgress
        (particle.phi(), secondCondBeta_, particle.omegaOU());
}
