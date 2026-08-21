/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SectionalParticleFlameSynthesis.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SectionalParticleFlameSynthesis<CloudType>::SectionalParticleFlameSynthesis
(
    const dictionary& dict,
    CloudType& owner
)
:
    SynthesisModel<CloudType>(dict,owner,typeName),
    
    speciesThermo_
    (
      dynamic_cast<const reactingMixture<ThermoType>&> (owner.thermo()).speciesData()
    ),
    
    inceptionRateIndex_(0),
    
    surfaceGrowthRateIndex_(1),
    
    primaryParticleNumberIndex_(2),
    
    particleNumberIndex_(3),
    
    geoMeanParticleVolumeIndex_(4),
    
    geoStdDevParticleVolumeIndex_(5),
    
    elsGasIndex_(6),
    
    elsParticleIndex_(7),
    
    elsIndex_(8),
        
    nPhysicalAerosolProperties_(9),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
    
    inception_(this->coeffDict().lookup("inception")),
    
    surfaceGrowth_(this->coeffDict().lookup("surfaceGrowth")),
    
    coagulation_(this->coeffDict().lookup("coagulation")),
    
    d00_(this->coeffDict().lookupOrDefault("d0", 1.0)),
    
    growth_(readScalar(this->coeffDict().lookup("growth"))),
    
    kf_(this->coeffDict().lookupOrDefault("kf", 1.0)),
    
    Df0_(readScalar(this->coeffDict().lookup("Df"))),
    
    coalescenceDiameter_(this->coeffDict().lookupOrDefault("coalescenceDiameter", 0.0)),
    
    dDfLimit_(this->coeffDict().lookupOrDefault("dLimit", VGREAT)),
    
    DfLimit_(this->coeffDict().lookupOrDefault("DfLimit", 3.0)),
    
    growthLimit_(this->coeffDict().lookupOrDefault("growthLimit", growth_)),

    elsGasNorm_(this->coeffDict().lookupOrDefault("ELSGasNorm", 1.0)),
    
    elsParticleNorm_(this->coeffDict().lookupOrDefault("ELSParticleNorm", 1.0)),

    zLower_(this->coeffDict().lookupOrDefault("zLower", 0.0)),

    zUpper_(this->coeffDict().lookupOrDefault("zUpper", 1.0)),
    
    d_(0.0),
    
    a_(0.0),
    
    v_(0.0),
    
    d0_(0.0),
    
    dc_(0.0),
    
    N_(0.0),
    
    Rg_(0.0),
    
    X_(0.0),
    
    beta_(0.0),

    precursorSpeciesDict_(this->coeffDict().subDict("precursorSpecies")),
    
    precursorSpeciesIndex_(-1),
    
    depositionSpeciesIndex_(-1),
    
    solidDensity_(readScalar(this->coeffDict().lookup("solidDensity"))),
    
    surfaceReactionEfficiency_(readScalar(this->coeffDict().lookup("surfaceReactionEfficiency"))),
    
    depositingSpeciesDict_(this->coeffDict().subDict("depositingSpecies")),
    
    coagulationModel_(this->coeffDict().lookup("coagulationModel")),
    
    turbCoagulation_(this->coeffDict().lookup("turbCoagulation")),

    coagulationSubSteps_(readLabel(this->coeffDict().lookup("coagulationSubSteps"))),

    mu_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("thermo:mu")),
    
    nut_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("nut")),
    
    epsilon_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("epsilon")),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    kB_(Ru_/NA_)
{
    Info << "Precursor species: " << endl;
    forAll (precursorSpeciesDict_.toc(), speciesI)
    {
        Info << token::TAB << precursorSpeciesDict_.toc()[speciesI] << endl;
    }    
    
    forAll(speciesThermo_, speciesI)
    {
        Info << speciesI << token::TAB << this->owner().slgThermo().carrier().species()[speciesI] << endl;
    }
    
    Info << "Solid density: "                   << solidDensity_        << endl;
    Info << "primary particle diameter: "       << d00_                 << endl;
    Info << "coalescence diameter: "            << coalescenceDiameter_ << endl;
    Info << "coalescence growth: "              << growth_              << endl;
    Info << "agglomeration growth: "            << growthLimit_         << endl;
    Info << "coalescence fractal dimesnion: "   << Df0_                 << endl;
    Info << "agglomeration fractal diemsnion: " << DfLimit_             << endl;
    Info << "Avogadro-Constant = "              << NA_                  << endl;
    Info << "Gas constant = "                   << Ru_                  << endl;
    Info << "Boltzman constant = "              << kB_                  << endl;
    
    Info << "Species for surface growth: " << endl;
    forAll (depositingSpeciesDict_.toc(), speciesI)
    {
        Info << token::TAB << depositingSpeciesDict_.toc()[speciesI] << endl;
    }
    
    Info << coagulationModel_ << " model is used as the coagulation model." << nl << endl;
    
    Info << coagulationSubSteps_ << " substeps applied for coagulation." << nl << endl;
    
    Info << "Normalisation factor for the gas based ELS signal: " << elsGasNorm_ << endl;
    Info << "Normalisation factor for the particle based ELS signal: " << elsParticleNorm_ << nl << endl;
    
    setSectionProperties();
    
    setSizeSplittingOperator();
}


template <class CloudType>
Foam::SectionalParticleFlameSynthesis<CloudType>::SectionalParticleFlameSynthesis
(
    const SectionalParticleFlameSynthesis<CloudType>& cm
)
:
    SynthesisModel<CloudType>(cm),
    
    speciesThermo_(cm.speciesThermo_),
        
    inceptionRateIndex_(0),
    
    surfaceGrowthRateIndex_(1),
    
    primaryParticleNumberIndex_(2),
    
    particleNumberIndex_(3),
    
    geoMeanParticleVolumeIndex_(4),
    
    geoStdDevParticleVolumeIndex_(5),
    
    elsGasIndex_(6),
    
    elsParticleIndex_(7),
    
    elsIndex_(8),
       
    nPhysicalAerosolProperties_(9),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
    
    inception_(this->coeffDict().lookup("inception")),
    
    surfaceGrowth_(this->coeffDict().lookup("surfaceGrowth")),
    
    coagulation_(this->coeffDict().lookup("coagulation")),
    
    d00_(this->coeffDict().lookupOrDefault("d0", 1.0)),
    
    growth_(readScalar(this->coeffDict().lookup("growth"))),
    
    kf_(this->coeffDict().lookupOrDefault("kf", 1.0)),
    
    Df0_(readScalar(this->coeffDict().lookup("Df"))),
    
    coalescenceDiameter_(this->coeffDict().lookupOrDefault("coalescenceDiameter", 0.0)),
    
    dDfLimit_(this->coeffDict().lookupOrDefault("dLimit", VGREAT)),
    
    DfLimit_(this->coeffDict().lookupOrDefault("DfLimit", 3.0)),
    
    growthLimit_(this->coeffDict().lookupOrDefault("growthLimit", growth_)),

    elsGasNorm_(this->coeffDict().lookupOrDefault("ELSGasNorm", 1.0)),
    
    elsParticleNorm_(this->coeffDict().lookupOrDefault("ELSParticleNorm", 1.0)),

    zLower_(this->coeffDict().lookupOrDefault("zLower", 0.0)),

    zUpper_(this->coeffDict().lookupOrDefault("zUpper", 1.0)),
    
    d_(0.0),
    
    a_(0.0),
    
    v_(0.0),
    
    d0_(0.0),
    
    dc_(0.0),
    
    N_(0.0),
    
    Rg_(0.0),
    
    X_(0.0),
    
    beta_(0.0),
            
    precursorSpeciesDict_(this->coeffDict().subDict("precursorSpecies")),
    
    solidSpeciesName_(this->coeffDict().lookup("solidSpeciesName")),

    precursorSpeciesIndex_(-1),
    
    depositionSpeciesIndex_(-1),
    
    solidDensity_(readScalar(this->coeffDict().lookup("solidDensity"))),
    
    surfaceReactionEfficiency_(readScalar(this->coeffDict().lookup("surfaceReactionEfficiency"))),
    
    depositingSpeciesDict_(this->coeffDict().subDict("depositingSpecies")),
    
    coagulationModel_(this->coeffDict().lookup("coagulationModel")),
    
    turbCoagulation_(this->coeffDict().lookup("turbCoagulation")),

    coagulationSubSteps_(readLabel(this->coeffDict().lookup("coagulationSubSteps"))),

    mu_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("thermo:mu")),
    
    nut_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("nut")),
    
    epsilon_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("epsilon")),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    kB_(Ru_/NA_)
{
    Info << "Precursor species: " << endl;
    forAll (precursorSpeciesDict_.toc(), speciesI)
    {
        Info << token::TAB << precursorSpeciesDict_.toc()[speciesI] << endl;
    }    
    
    Info << "Solid density: "                   << solidDensity_        << endl;
    Info << "primary particle diameter: "       << d00_                 << endl;
    Info << "coalescence diameter: "            << coalescenceDiameter_ << endl;
    Info << "coalescence growth: "              << growth_              << endl;
    Info << "agglomeration growth: "            << growthLimit_         << endl;
    Info << "coalescence fractal dimesnion: "   << Df0_                 << endl;
    Info << "agglomeration fractal diemsnion: " << DfLimit_             << endl;
    Info << "diameter growth: "                 << growth_              << endl;
    Info << "Avogadro-Constant = "              << NA_                  << endl;
    Info << "Gas constant = "                   << Ru_                  << endl;
    Info << "Boltzman constant = "              << kB_                  << endl;
    
    Info << "Species for surface growth: " << endl;
    forAll (depositingSpeciesDict_.toc(), speciesI)
    {
        Info << token::TAB << depositingSpeciesDict_.toc()[speciesI] << endl;
    }
    
    Info << coagulationModel_ << " model is used as the coagulation model." << nl << endl;

    Info << coagulationSubSteps_ << " substeps applied for coagulation." << nl << endl;
    
    Info << "Normalisation factor for the gas based ELS signal: " << elsGasNorm_ << endl;
    Info << "Normalisation factor for the particle based ELS signal: " << elsParticleNorm_ << nl << endl;
        
    setSectionProperties();
    
    setSizeSplittingOperator();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SectionalParticleFlameSynthesis<CloudType>::~SectionalParticleFlameSynthesis()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::synthesize()
{
    if (inception_)
    {
        Info << "Solid particle inception..." << endl;    
    
        forAllIters(this->owner(), particle)
        {
            particleInception(particle());
        }
    }
    
    if (surfaceGrowth_)
    {
        Info << "Solid particle surface growth..." << endl;    
        
        forAllIters(this->owner(), particle)
        {
            particleSurfaceGrowth(particle());
        }
    }
    
    if (coagulation_)
    {
        Info << "Solid particle coagulation..." << endl;
        
        forAllIters(this->owner(), particle)
        {
            particleCoagulation(particle());
        }
    }
    
    forAllIters(this->owner(), particle)
    {
        calcPhysicalAerosolProperties(particle());
    }
    
    Info << "Finished!" << endl;    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::setSectionProperties()
{   
    const label nSecSqr = nSections_*nSections_;

    d_.setSize (nSections_+1, 0.0);
    a_.setSize (nSections_+1, 0.0);
    v_.setSize (nSections_+1, 0.0);
    d0_.setSize(nSections_+1, 0.0);
    dc_.setSize(nSections_+1, 0.0);
    N_.setSize (nSections_+1, 0.0);
    Df_.setSize(nSections_+1, 0.0);
    Rg_.setSize(nSections_+1, 0.0);
    
    v_[0] = M_PI/6.0*Foam::pow(d00_, 3.0);
    
    const scalar& thetaDeg = 90;
    const scalar& thetaRad = thetaDeg*M_PI/180;
    const scalar& q = 4.0*M_PI/532e-09*Foam::sin(thetaRad/2.0);
    
    Info << "Volumetric spacing factor is " << growth_ << nl << endl;
    
    Info << "section"          << token::TAB 
         << "diameter"         << token::TAB 
         << "primary-diameter" << token::TAB 
         << "surface-area"     << token::TAB 
	 << "volume"           << token::TAB 
         << "N"                << token::TAB 
         << "Df"               << token::TAB 
         << "dc"               << token::TAB 
         << "Rg"               << token::TAB 
         << "qRg"              << endl;

    scalar Ntmp0 = 0.0;
    label transition = 0;
    scalar growth = growthLimit_;

    forAll (v_, section)
    {
        scalar Ntmp = floor(Ntmp0*growth);
        if(Ntmp == Ntmp0)
            Ntmp += 1;     
	v_[section] = v_[transition]*Ntmp;
        d_[section] = Foam::pow((6.0*v_[section]/M_PI),1.0/3.0);
        a_[section] = M_PI*Foam::pow(d_[section], 2.0);

	if(d_[section]>coalescenceDiameter_ && transition==0)
        {
            Ntmp = ceil(M_PI/6.0*Foam::pow(coalescenceDiameter_,3.0)/v_[0]);
            v_[section] = Ntmp*v_[0];
            d_[section] = Foam::pow((6.0*v_[section]/M_PI),1.0/3.0);
            a_[section] = M_PI*Foam::pow(d_[section], 2.0);

            Ntmp=1.0;
            transition=section;
        }

        Df_[section]=DfLimit_;
        if(Ntmp==1.0)
            Df_[section] = 3.0;
        else if(d_[section]>coalescenceDiameter_)
        {
            Df_[section] = Df0_;
            growth = growth_;
        }

        d0_[section] = d_[transition];
        N_[section] = Ntmp;
        Ntmp0 = Ntmp;
           
        dc_[section] = d0_[section]*Foam::pow(N_[section], 1.0/Df_[section]);    
        Rg_[section] = Foam::pow(N_[section]/kf_, 1.0/Df_[section])*(d0_[section]/2.0);
        
        Info << section        << token::TAB 
             << d_[section]    << token::TAB 
             << d0_[section]   << token::TAB 
             << a_[section]    << token::TAB 
             << v_[section]    << token::TAB 
             << N_[section]    << token::TAB 
             << Df_[section]   << token::TAB 
             << dc_[section]   << token::TAB 
             << Rg_[section]   << token::TAB 
             << q*Rg_[section] << endl;
    }
    
    beta_.setSize(nSecSqr, 0.0);
}


template<class CloudType>
const Foam::wordList Foam::SectionalParticleFlameSynthesis<CloudType>::psdPropertyNames() const
{
    wordList sectionNames;
    sectionNames.setSize(nSections_);
    
    forAll(sectionNames,section)
    {
        word sectionName;
        
        if (section < 10)
        {
            sectionNames[section] = "section" + Foam::name(0) + Foam::name(section);
        }
        else
        {
            sectionNames[section] = "section" + Foam::name(section);
        }
    }
        
    return sectionNames;
}


template<class CloudType>
const Foam::wordList Foam::SectionalParticleFlameSynthesis<CloudType>::physicalAerosolPropertyNames() const
{
    wordList physicalAerosolPropertyNames;
    physicalAerosolPropertyNames.setSize(nPhysicalAerosolProperties_);
    
    physicalAerosolPropertyNames[0]  = "inceptionRate";
    physicalAerosolPropertyNames[1]  = "surfaceGrowthRate";
    physicalAerosolPropertyNames[2]  = "primaryParticleNumber";
    physicalAerosolPropertyNames[3]  = "particleNumber";
    physicalAerosolPropertyNames[4]  = "geoMeanParticleVolume";
    physicalAerosolPropertyNames[5]  = "geoStdDevParticleVolume";
    physicalAerosolPropertyNames[6]  = "ELSgas";
    physicalAerosolPropertyNames[7]  = "ELSparticle";
    physicalAerosolPropertyNames[8]  = "ELS";
            
    return physicalAerosolPropertyNames;
}


template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::setSizeSplittingOperator()
{
    const label nSecSqr = nSections_*nSections_;
    const label nSecPow3 = nSecSqr*nSections_;

    X_.setSize(nSecPow3, 0.0);
    
    for (label i=0; i <= nSections_-1; i++)
    {
        for (label j=0; j <= nSections_-1; j++)
        {
            for (label k=1; k <= nSections_-1; k++)
            {
                // conditions in parentheses check if the combined volume of colliding particles i and j is between k and k+1.
                if (v_[k] <= (v_[i]+v_[j]) && (v_[i]+v_[j]) < v_[k+1])
                {
                    X_[i*nSecSqr+j*nSections_+k] = (v_[k+1]-v_[i]-v_[j])/(v_[k+1]-v_[k]);
                }             
                else
                {
                    if (v_[k-1] <= (v_[i]+v_[j]) && (v_[i]+v_[j]) < v_[k])
                    {
                        X_[i*nSecSqr+j*nSections_+k] = (v_[i]+v_[j]-v_[k-1])/(v_[k]-v_[k-1]);
                    }
                    else 
                    {
                        X_[i*nSecSqr+j*nSections_+k] = 0;
                    }
                }              
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::particleInception(particleType& particle) 
{
    CloudType& cloud(this->owner());
    const scalar deltaT = cloud.mesh().time().deltaT().value(); // [s]

    const scalar molecularWeightMixture = calcMolecularWeightMixture(particle.Y());   // [kg/kmol]

    particle.physicalAerosolProperties()[inceptionRateIndex_] = 0.0;
    
    forAll(speciesThermo_, speciesI)
    {
        const word speciesName = this->owner().slgThermo().carrier().species()[speciesI];
 
        if (precursorSpeciesDict_.found(speciesName))
        {
            precursorSpeciesIndex_ = speciesI;
 
            //const scalar molecularWeightPrecursorSpecies = speciesThermo_[precursorSpeciesIndex_].W(); // [kg/kmol]
 
            // calculate the number of particles inserted into first section du to formation of precursor species
            //const scalar insertedParticles = NA_*particle.Y()[precursorSpeciesIndex_]*particle.pc()*molecularWeightMixture/(Ru_*particle.T()*molecularWeightPrecursorSpecies); // [#/m^3] = [1/kmol] [1] [kg/m s^2] [kg/kmol] / [kg m^2/s^2 kmol K] [K] [kg/kmol]
            const scalar insertedParticles = 6.0*particle.Y()[precursorSpeciesIndex_]*particle.pc()*molecularWeightMixture / (M_PI*Foam::pow(d00_,3.0)*solidDensity_*Ru_*particle.T());
             
            particle.psdProperties()[0] += insertedParticles;
         
            // calculate the nucleation rate
            particle.physicalAerosolProperties()[inceptionRateIndex_] += insertedParticles/deltaT; // [#/m^3 s] 

            scalar dPrecursorSpecie = -particle.Y()[precursorSpeciesIndex_];
            scalar dm = dPrecursorSpecie * particle.m();

            particle.dpYsource()[precursorSpeciesIndex_] = dPrecursorSpecie;
            particle.dpMsource() = dm;
            //TODO: Eulerian mass field should be updated too (has small impact thou)
            particle.updateDispersedSources();
            particle.wt()=particle.m()/cloud.deltaM();
        }
    }
}


template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::particleSurfaceGrowth(particleType& particle) 
{
    if ((particle.XiC()[0] <= zLower_) || (particle.XiC()[0] >= zUpper_))
       return;

    CloudType& cloud(this->owner());
    const scalar deltaT = cloud.mesh().time().deltaT().value(); // [s]

    const scalar molecularWeightMixture = calcMolecularWeightMixture(particle.Y()); // [kg/kmol]
    
    const scalar densityMixture = particle.pc()*molecularWeightMixture/(Ru_*particle.T()); // [kg/m^3] = [kg/m s^2] [kg/kmol] / [kg m^2/s^2 kmol K] [K]
    
    particle.physicalAerosolProperties()[surfaceGrowthRateIndex_] = 0.0;
        
    forAll(speciesThermo_, speciesI)
    { 
        const word speciesName = this->owner().slgThermo().carrier().species()[speciesI];
        
        if (depositingSpeciesDict_.found(speciesName))
        {            
            depositionSpeciesIndex_ = speciesI;
            
            scalarList sectionChange(nSections_, 0.0);

            const scalar molecularWeightDepositingSpecies = speciesThermo_[depositionSpeciesIndex_].W();                           // [kg/kmol]

            const scalar monomerMassDepositingSpecies     = molecularWeightDepositingSpecies/NA_;                                  // [kg/#] = [kg/kmol] / [#/kmol]
            const scalar monomerVolumeDepositingSpecies   = monomerMassDepositingSpecies/solidDensity_;                            // [m^3/#] = [kg/#] / [kg/m^3]
            
            forAll (particle.psdProperties(), section)
            {
                // calculation of the volumetric growth via the free molecular collision kernel approach
                const scalar d                     = Foam::pow(6.0*monomerVolumeDepositingSpecies/M_PI,1.0/3.0);                   // [m]
                   
                const scalar gasPhaseConcentration = particle.Y()[depositionSpeciesIndex_]*densityMixture/molecularWeightDepositingSpecies; // [kmol/m^3]
                
                const scalar Rsurfrxn              = 2.2*surfaceReactionEfficiency_*gasPhaseConcentration*Foam::sqrt((M_PI*kB_*particle.T())/(2*monomerMassDepositingSpecies))*Foam::pow(d + dc_[section], 2.0); // [kmol/s] per particle
                
                const scalar dmfm                  = Rsurfrxn*molecularWeightDepositingSpecies*deltaT;                             // [kg] = [kmol/s] [kg/kmol] [s] per particle

                // calculation of the volumetric growth via the continuum collision kernel approach
                const scalar sigma                 = 0.5*(3.667e-10 + d);
                const scalar eps                   = Foam::sqrt(190.0*97); // 190 for CO2 due to the similarity to SiO2, 97 due to the strong presence of N2
                const scalar keps                  = kB_/eps; 
                const scalar Tstar                 = particle.T()*keps; 
                const scalar omega                 = 1.06036/Foam::pow(Tstar, 0.15610) + 0.19300/Foam::exp(0.47635*Tstar) + 1.03587/Foam::exp(1.52995*Tstar) + 1.76474/Foam::exp(3.89411*Tstar);
                const scalar D                     = 1.883e-26*Foam::pow(particle.T(), 1.5)/(particle.pc()*Foam::pow(sigma, 2.0)*omega)*Foam::sqrt(1.0/molecularWeightDepositingSpecies + 1.0/molecularWeightMixture); // [m^2/s]
                const scalar aggSurf               = a_[0]*(section+1);
                
                const scalar dmco                  = 2*D*particle.Y()[depositionSpeciesIndex_]/d_[section]*aggSurf*solidDensity_*deltaT; // [kg] = [m^2/s] [-] [m^2] [s] [kg/m^3] / [m]

                // calculation of the volumetric growth due to the harmonic mean
                scalar dm                          = (dmfm*dmco)/(dmfm + dmco + SMALL);                                            // [kg]
                
                scalar depositedMassFraction       = dm*particle.psdProperties()[section]/densityMixture;                          // [-] = [kg] [#/m^3] / [kg/m^3]
                
                if (depositedMassFraction > particle.Y()[depositionSpeciesIndex_])
                {
                    dm *= particle.Y()[depositionSpeciesIndex_]/depositedMassFraction;
                    depositedMassFraction = particle.Y()[depositionSpeciesIndex_];
                }
                
                const scalar dV = dm/solidDensity_;                                                                                // [m^3] = [-] [kg] / [kg/m^3]
                
                particle.physicalAerosolProperties()[surfaceGrowthRateIndex_] += dV/deltaT;
              
                // distribution of number density        
                scalar targetVolume = v_[section] + dV;                                                                            // [m^3] = [m^3] + [m^3]
 
                if (targetVolume > v_[nSections_-1])
                {
                    targetVolume = v_[nSections_-1];
                }
                
                label targetSectionL = nSections_-2;
                label targetSectionU = nSections_-1;
                    
                for (label searchSection = 0; searchSection < nSections_-2; searchSection++)
                {
                    if ((targetVolume >= v_[searchSection]) && (targetVolume <= v_[searchSection+1]))
                    {
                        targetSectionL = searchSection;
                        targetSectionU = searchSection + 1;
                    }
                }            

                const scalar xiU = (targetVolume-v_[targetSectionL])/(v_[targetSectionU] - v_[targetSectionL]);
                const scalar xiL = 1 - xiU;
                    
                const scalar nParticlesToEnterL = xiL*particle.psdProperties()[section];
                const scalar nParticlesToEnterU = xiU*particle.psdProperties()[section];
                const scalar nParticlesToLeave  = nParticlesToEnterL + nParticlesToEnterU;
                    
                sectionChange[targetSectionL] += nParticlesToEnterL;
                sectionChange[targetSectionU] += nParticlesToEnterU;
                sectionChange[section]        -= nParticlesToLeave;

                particle.Y()[depositionSpeciesIndex_]    = max(particle.Y()[depositionSpeciesIndex_] - depositedMassFraction, 0.0); 
            }
            
            forAll (particle.psdProperties(), section)
            {
                particle.psdProperties()[section] += sectionChange[section];
            }
        }
    }
    
    particle.physicalAerosolProperties()[surfaceGrowthRateIndex_] /= nSections_;
}

template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::particleCoagulation(particleType& particle) 
{
    scalar particleNumber = 0.0;

    forAll (particle.psdProperties(), section)
        particleNumber += particle.psdProperties()[section];

    if ((particle.XiC()[0] < zLower_) || (particle.XiC()[0] > zUpper_))// || particleNumber < 2.0)
       return;

    CloudType& cloud(this->owner());
    scalar deltaT = cloud.mesh().time().deltaT().value();
    
    const label nSecSqr = nSections_*nSections_;
   
    if (coagulationModel_ == "freeMolecularCollisionKernel")
    {
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                const scalar v1 =       1.0/v_[i]           + 1.0/v_[j];
                const scalar v2 = Foam::pow(v_[i], 1.0/3.0) + pow(v_[j], 1.0/3.0);
                
                beta_[i*nSections_+j] = Foam::pow(3.0/(4.0*M_PI), 1.0/6.0)*Foam::pow((6.0*kB_*particle.T()/solidDensity_), 0.5)*pow(v1, 0.5)*pow(v2, 2.0);               
            }
        }
    } 
    else if (coagulationModel_ == "freeMolecularCollisionKernelFriedlander")
    {
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                const scalar v1 = 1.0/v_[i] + 1.0/v_[j];
                const scalar v2 = Foam::pow(v_[i], 1.0/Df_[i]) + pow(v_[j], 1.0/Df_[j]);
                
                const scalar lambda = 2.0/Df0_ - 0.5;
                const scalar C = Foam::pow(dc_[j]/2.0, 2-6.0/Df0_); // caution dc_[j]
                
                beta_[i*nSections_+j] = Foam::pow((6.0*kB_*particle.T()/solidDensity_), 0.5)*Foam::pow(3.0/(4.0*M_PI), lambda)*C*pow(v1, 0.5)*pow(v2, 2.0);               
            }
        }
    }
    else if (coagulationModel_ == "freeMolecularCollisionKernelDf")
    {
        interpolationCellPoint<scalar> nutIntp_(this->nut_);
        const scalar& nut = nutIntp_.interpolate(particle.position(),particle.cell(),particle.face());

        interpolationCellPoint<scalar> epsIntp_(this->epsilon_);
        const scalar& eps = epsIntp_.interpolate(particle.position(),particle.cell(),particle.face());

        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                const scalar v1 =       1.0/N_[i]            +       1.0/N_[j];
                const scalar v2 = Foam::pow(N_[i], 1.0/Df0_) + Foam::pow(N_[j], 1.0/Df0_);

                const scalar C = 4.89;
                const scalar a  = Foam::pow((3.0*v_[0])/(4.0*M_PI), 1.0/6.0)*Foam::pow((6.0*kB_*particle.T()/solidDensity_), 0.5); // caution !!! v[0]
                const scalar a1 = Foam::pow(2.0, Df0_)*a/C;

                scalar beta = a1*Foam::pow(v1, 0.5)*Foam::pow(v2, Df0_);

                if(turbCoagulation_)
                    beta = max(beta,1.3*Foam::pow(eps/nut, 0.5)*Foam::pow(Foam::pow(d_[i]/2.0, 1.0/Df0_) + Foam::pow(d_[j]/2.0, 1.0/Df0_), Df0_));

                beta_[i*nSections_+j] = beta;
            }
        }
    }
    else if (coagulationModel_ == "continuum")
    {
        scalar A = 1.966e-06;
        scalar B = 147.47;
        
        const scalar mu = A*Foam::pow(particle.T(), 1.5)/(particle.T() + B);
        
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                // continuum regime
                const scalar v1 = 1.0/Foam::pow(v_[i], 1.0/Df_[i]) + 1.0/Foam::pow(v_[j], 1.0/Df_[j]);
                const scalar v2 =     Foam::pow(v_[i], 1.0/Df_[i]) +     Foam::pow(v_[j], 1.0/Df_[j]);
                                
                beta_[i*nSections_+j] = (2.0*kB_*particle.T())/(3.0*mu)*v1*v2;
            }
        }
    }
    else if (coagulationModel_ == "FuchsSutuginCollisionKernel")
    {
        interpolationCellPoint<scalar> muIntp_(this->mu_);
        const scalar& mu = muIntp_.interpolate(particle.position(),particle.cell(),particle.face());   
       
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                const scalar dci = dc_[i];
                const scalar dcj = dc_[j];
                
                const scalar mp1 = v_[i]*solidDensity_;
                const scalar mp2 = v_[j]*solidDensity_;
      
                const scalar lambda1 = mu/particle.pc()*sqrt((M_PI*kB_*particle.T())/(2*mp1));
                const scalar lambda2 = mu/particle.pc()*sqrt((M_PI*kB_*particle.T())/(2*mp2));
                
                const scalar Kn1 = (2.0*lambda1)/dci;
                const scalar Kn2 = (2.0*lambda2)/dcj;
                
                const scalar D1 = (kB_*particle.T())/(3.0*M_PI*mu*dci)*((5.0 + 4.0*Kn1 + 6.0*Kn1*Kn1 + 18.0*Kn1*Kn1*Kn1)/(5.0 - Kn1 + (8.0 + M_PI)*Kn1*Kn1));
                const scalar D2 = (kB_*particle.T())/(3.0*M_PI*mu*dcj)*((5.0 + 4.0*Kn2 + 6.0*Kn2*Kn2 + 18.0*Kn2*Kn2*Kn2)/(5.0 - Kn2 + (8.0 + M_PI)*Kn2*Kn2));
                
                const scalar c1 = sqrt((8.0*kB_*particle.T())/(M_PI*mp1));
                const scalar c2 = sqrt((8.0*kB_*particle.T())/(M_PI*mp2));
                
                const scalar l1 = (8.0*D1)/(M_PI*c1);
                const scalar l2 = (8.0*D2)/(M_PI*c2);
                
                const scalar g1 = (Foam::pow((dci + l1), 3) - Foam::pow((dci*dci + l1*l1), 1.5))/(3.0*dci*l1) - dci;
                const scalar g2 = (Foam::pow((dcj + l2), 3) - Foam::pow((dcj*dcj + l2*l2), 1.5))/(3.0*dcj*l2) - dcj;

                scalar beta = 2.0*M_PI*(D1 + D2)*(dci + dcj)/((dci + dcj)/(dci + dcj + 2.0*sqrt(g1*g1 + g2*g2)) + (8.0*(D1 + D2))/(sqrt(c1*c1 + c2*c2)*(dci + dcj)));

                beta_[i*nSections_+j] = beta;
            }
        }
    }
    else if (coagulationModel_ == "FuchsSutuginCollisionKernel"  && turbCoagulation_)
    {
        interpolationCellPoint<scalar> muIntp_(this->mu_);
        const scalar& mu = muIntp_.interpolate(particle.position(),particle.cell(),particle.face());   

        interpolationCellPoint<scalar> nutIntp_(this->nut_);
        const scalar& nut = nutIntp_.interpolate(particle.position(),particle.cell(),particle.face());

        interpolationCellPoint<scalar> epsIntp_(this->epsilon_);
        const scalar& eps = epsIntp_.interpolate(particle.position(),particle.cell(),particle.face());
        
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                const scalar dci = dc_[i];
                const scalar dcj = dc_[j];
                
                const scalar mp1 = v_[i]*solidDensity_;
                const scalar mp2 = v_[j]*solidDensity_;
      
                const scalar lambda1 = mu/particle.pc()*sqrt((M_PI*kB_*particle.T())/(2*mp1));
                const scalar lambda2 = mu/particle.pc()*sqrt((M_PI*kB_*particle.T())/(2*mp2));

                const scalar Kn1 = (2.0*lambda1)/dci;
                const scalar Kn2 = (2.0*lambda2)/dcj;

                const scalar D1 = (kB_*particle.T())/(3.0*M_PI*mu*dci)*((5.0 + 4.0*Kn1 + 6.0*Kn1*Kn1 + 18.0*Kn1*Kn1*Kn1)/(5.0 - Kn1 + (8.0 + M_PI)*Kn1*Kn1));
                const scalar D2 = (kB_*particle.T())/(3.0*M_PI*mu*dcj)*((5.0 + 4.0*Kn2 + 6.0*Kn2*Kn2 + 18.0*Kn2*Kn2*Kn2)/(5.0 - Kn2 + (8.0 + M_PI)*Kn2*Kn2));

                const scalar c1 = sqrt((8.0*kB_*particle.T())/(M_PI*mp1));
                const scalar c2 = sqrt((8.0*kB_*particle.T())/(M_PI*mp2));

                const scalar l1 = (8.0*D1)/(M_PI*c1);
                const scalar l2 = (8.0*D2)/(M_PI*c2);

                const scalar g1 = (Foam::pow((dci + l1), 3) - Foam::pow((dci*dci + l1*l1), 1.5))/(3.0*dci*l1);// - dci;
                const scalar g2 = (Foam::pow((dcj + l2), 3) - Foam::pow((dcj*dcj + l2*l2), 1.5))/(3.0*dcj*l2);// - dcj;

                const scalar betaLami = 2.0*M_PI*(D1 + D2)*(dci + dcj)/((dci + dcj)/(dci + dcj + 2.0*sqrt(g1*g1 + g2*g2)) + (8.0*(D1 + D2))/(sqrt(c1*c1 + c2*c2)*(dci + dcj)));

                const scalar betaTurb = 1.3*Foam::pow(eps/nut, 0.5)*Foam::pow(Foam::pow(dci/2.0, 1.0/Df_[i]) + Foam::pow(dcj/2.0, 1.0/Df_[j]), (Df_[i]+Df_[j])/2.0);

                beta_[i*nSections_+j] = max(betaLami,betaTurb);

                if (betaTurb > betaLami)
                    Info << "betaLami = " << betaLami << token::TAB << "betaTurb = " << betaTurb << endl;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::SectionalParticleFlameSynthesis<CloudType>::particleCoagulation(particleType& particle)"
        )   << "no coagulation model set" 
            << exit(FatalError);
    }

//// Test number conserved //////////////////////////////////////
    scalar totalVolumeBefore = 0.0;//////////////////////////////
    for (label k=0; k <= nSections_-1; k++)//////////////////////
        totalVolumeBefore += v_[k]*particle.psdProperties()[k];//
/////////////////////////////////////////////////////////////////

    scalar subDeltaT = deltaT/coagulationSubSteps_;
    
    for (label subStep=1; subStep <= coagulationSubSteps_; subStep++)
    {    
        scalarField psdProperties_tmp = particle.psdProperties();

        // calculating the gain and loss terms due to coagulation.
        for (label k=0; k <= nSections_-1; k++)
        {
            scalar add = 0.0; // addition term when i and j collide to form a k sized particle.
            scalar sub = 0.0; // subtraction term when k collides with any other particle.

            for (label i=0; i <= nSections_-1; i++)
            {
                sub = sub + beta_[k*nSections_+i]*particle.psdProperties()[i];

                if(i<=k)
                {
                    for (label j=0; j<=i; j++)
                    {
                        scalar X = X_[i*nSecSqr+j*nSections_+k];
                        if(X!=0.0)
                        {
                            scalar notTwice = 1.0;
                            if(i==j) 
                                notTwice =.5;
                            add = add + (X*beta_[i*nSections_+j]*particle.psdProperties()[i]*particle.psdProperties()[j]*notTwice);
                        }
                    }
                }
            }

            sub*=particle.psdProperties()[k];

            psdProperties_tmp[k] += subDeltaT*(add - sub);
	}
        particle.psdProperties() = psdProperties_tmp;
    }
//// Test number conserved /////////////////////////////////////////////////////////////
    scalar totalVolumeAfter = 0.0;//////////////////////////////////////////////////////
    for (label k=0; k <= nSections_-1; k++)/////////////////////////////////////////////
        totalVolumeAfter += v_[k]*particle.psdProperties()[k];//////////////////////////
    if(totalVolumeBefore>0.0)///////////////////////////////////////////////////////////
        if(mag(totalVolumeAfter-totalVolumeBefore)/totalVolumeBefore>0.001)/////////////
            Pout << "Number conergence ERROR! n-before: " << totalVolumeBefore//////////
                 << "   n-after: "      << totalVolumeAfter/////////////////////////////
                 << "   n-difference: " << totalVolumeAfter-totalVolumeBefore << endl;//
////////////////////////////////////////////////////////////////////////////////////////
}


template <class CloudType>
Foam::scalar Foam::SectionalParticleFlameSynthesis<CloudType>::calcMolecularWeightMixture(scalarField& Y) 
{
    scalar molecularWeightMixture = 0.0;
    
    forAll(speciesThermo_,speciesI) 
    {
        molecularWeightMixture += Y[speciesI]/speciesThermo_[speciesI].W();
    }
        
    return 1.0/molecularWeightMixture;    
}


template <class CloudType>
Foam::scalarList Foam::SectionalParticleFlameSynthesis<CloudType>::X(scalarField& Y)
{
    scalar molecularWeightMixture = 0.0;

    forAll(speciesThermo_,speciesI)
    {
        molecularWeightMixture += Y[speciesI]/speciesThermo_[speciesI].W();
    }

    molecularWeightMixture = 1.0/molecularWeightMixture;

    scalarList X(Y.size(), 0.0);

    forAll(speciesThermo_,speciesI)
    {
        X[speciesI] = Y[speciesI]*speciesThermo_[speciesI].W()/molecularWeightMixture;
    }

    return X;
}


template <class CloudType>
void Foam::SectionalParticleFlameSynthesis<CloudType>::calcPhysicalAerosolProperties(particleType& particle) 
{
    // elastic light scattering signal of the gas phase
    scalarList x = X(particle.Y());

    scalar CSM = 0.0;
    forAll(speciesThermo_, speciesI)
    {
        const word speciesName = this->owner().slgThermo().carrier().species()[speciesI];
        if (speciesName == "CO2")
            CSM+=2.3907*x[speciesI];
	else if (speciesName == "O2")
            CSM+=0.8592*x[speciesI];
        else if (speciesName == "CO")
            CSM+=1.2446*x[speciesI];
        else if (speciesName == "N2")
            CSM+=1.0000*x[speciesI];
        else if (speciesName == "CH4")
            CSM+=2.1337*x[speciesI];
        else if (speciesName == "H2O")
            CSM+=0.6946*x[speciesI];
        else if (speciesName == "H2")
            CSM+=0.2124*x[speciesI];
        else if (speciesName == "OH")
            CSM+=1.4859*x[speciesI];
        else if (speciesName == "Ar")
            CSM+=0.8650*x[speciesI];
        else if (speciesName == "O")
            CSM+=0.1713*x[speciesI];
        else if (speciesName == "H")
            CSM+=0.1479*x[speciesI];
        else if (speciesName == "C2H6O")
            CSM+=8.6473*x[speciesI];
        else if (speciesName == "C2H4")
            CSM+=5.8029*x[speciesI];
        else if (speciesName == "C2H6")
            CSM+=6.2984*x[speciesI];
        else if (speciesName == "C3H8")
            CSM+=12.7542*x[speciesI];
        else if (speciesName == "CH2O")
            CSM+=0.0*x[speciesI];
        else if (speciesName == "C2H2")
            CSM+=4.0096*x[speciesI];
        else if (speciesName == "CH3")
            CSM+=1.5770*x[speciesI];
        else if (speciesName == "He")
            CSM+=0.0132*x[speciesI];
        else if (speciesName == "NO")
            CSM+=0.9834*x[speciesI];
        else if (speciesName == "C2H5OH")
            CSM+=8.0026*x[speciesI];
        else if (speciesName == "C2H4O")
            CSM+=6.7971*x[speciesI];
        else if (speciesName == "C3H6O")
            CSM+=13.2247*x[speciesI];
    }

    particle.physicalAerosolProperties()[elsGasIndex_] = CSM/particle.T();
   
    // elastic light scattering signal of the particle phase
    particle.physicalAerosolProperties()[elsParticleIndex_]             = 0.0;
    
    const scalar lambda = 532e-09;
    const scalar deg = 90;
    const scalar rad = deg*M_PI/180;
    const scalar q = 4.0*M_PI/lambda*Foam::sin(rad/2.0);
    const scalar k = 2*M_PI/lambda;
    
    forAll (particle.psdProperties(), section)
    {
        scalar Caggsca = 0.0;
                
        const scalar qR = q*Rg_[section];

        if (qR < 0.1)
            //Caggsca = Foam::pow(d_[0], 6 - 2*Df_[section])*Foam::pow(Rg_[section], 2*Df_[section]);
            Caggsca = Foam::pow(N_[section], 2)*Foam::pow(k, 4)*Foam::pow(d0_[section]/2, 6);
        else if ((qR > 0.1) && (qR < 1))
            //Caggsca = Foam::pow(d_[0], 6 - 2*Df_[section])*(Foam::pow(Rg_[section], 2*Df_[section]) - Foam::pow(Rg_[section], 4*Df_[section]));
            Caggsca = Foam::pow(N_[section], 2)*Foam::pow(k, 4)*Foam::pow(d0_[section]/2, 6)*(1-(Foam::pow(q, 2)/3*Foam::pow(Rg_[section], 2)));
        else if (qR > 1)
            //Caggsca = Foam::pow(d_[0], 6 - 2*Df_[section])*Foam::pow(Rg_[section], Df_[section]);
            Caggsca = Foam::pow(N_[section], 2)*Foam::pow(k, 4)*Foam::pow(d0_[section]/2, 6)*0.77*Foam::pow(qR, -Df_[section]);
        
        particle.physicalAerosolProperties()[elsParticleIndex_]           += particle.psdProperties()[section]*Caggsca;
    }
    
    particle.physicalAerosolProperties()[elsIndex_] = particle.physicalAerosolProperties()[elsGasIndex_]/elsGasNorm_ + particle.physicalAerosolProperties()[elsParticleIndex_]/elsParticleNorm_;
    
    // calculation of the particle number density
    particle.physicalAerosolProperties()[particleNumberIndex_] = 0.0;
    
    forAll (particle.psdProperties(), section)
    {
        particle.physicalAerosolProperties()[particleNumberIndex_] += particle.psdProperties()[section];
    }
    
    // calculation of the primary particle number density
    particle.physicalAerosolProperties()[primaryParticleNumberIndex_] = 0.0;
    
    forAll (particle.psdProperties(), section)
    {
        particle.physicalAerosolProperties()[primaryParticleNumberIndex_] += particle.psdProperties()[section]*N_[section];
    }
}

// ************************************************************************* //
