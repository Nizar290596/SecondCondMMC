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

#include "SectionalHomogeneousNucleationCondensation.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SectionalHomogeneousNucleationCondensation<CloudType>::SectionalHomogeneousNucleationCondensation
(
    const dictionary& dict,
    CloudType& owner
)
:
    SynthesisModel<CloudType>(dict,owner,typeName),
    
    specieThermo_
    (
      dynamic_cast<const reactingMixture<ThermoType>&> (owner.thermo()).speciesData()
    ),
    
    saturationIndex_(0),
    
    nucleationRateIndex_(1),
    
    condensationRateIndex_(2),
    
    nPhysicalAerosolProperties_(3),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
    
    nuclSpecies_(this->coeffDict().lookup("nucleatingSpecies")),

    inertSpecies_(this->coeffDict().lookup("inertSpecies")),
    
    surfaceTensionModelCoeffs_(this->coeffDict().subDict("surfaceTensionModelCoeffs")),

    surfaceTensionModelCoeffM_(readScalar(surfaceTensionModelCoeffs_.lookup("m"))),

    surfaceTensionModelCoeffN_(readScalar(surfaceTensionModelCoeffs_.lookup("n"))),

    liquidDensityModelCoeffs_(this->coeffDict().subDict("liquidDensityModelCoeffs")),

    liquidDensityModelCoeffM_(readScalar(liquidDensityModelCoeffs_.lookup("m"))),

    liquidDensityModelCoeffN_(readScalar(liquidDensityModelCoeffs_.lookup("n"))),

    saturationPressureModelCoeffs_(this->coeffDict().subDict("saturationPressureModelCoeffs")),

    saturationPressureModelCoeffA_(readScalar(saturationPressureModelCoeffs_.lookup("a"))),

    saturationPressureModelCoeffB_(readScalar(saturationPressureModelCoeffs_.lookup("b"))),

    saturationPressureModelCoeffC_(readScalar(saturationPressureModelCoeffs_.lookup("c"))),   
 
    nucleationRateModel_(this->coeffDict().lookup("nucleationRateModel")),
    
    condensationRateModel_(this->coeffDict().lookup("condensationRateModel")),
    
    d_(),
    
    a_(),
    
    v_(),
    
    nuclSpeciesIndex_(0),
    
    inertSpeciesIndex_(0),
    
    surfaceTension_(0),
    
    liquidDensity_(0),
    
    molecularWeightNuclSpecies_(0),
    
    molecularMass_(0),
    
    molecularVolume_(0),
    
    molecularSurface_(0),
    
    molecularWeightMixture_(0),
    
    saturationPressure_(0),
    
    vaporPressure_(0),
    
    mixtureDensity_(0),
    
    nucleationRate_(0),
    
    condensationRate_(0),
    
    meanFreePath_(0),
    
    d0_(readScalar(this->coeffDict().lookup("d0"))),
       
    volGrowthFactor_(readScalar(this->coeffDict().lookup("volGrowthFactor"))),
    
    criticalDiameter_(0),
    
    criticalVolume_(0),
    
    condensedMassFraction_(0),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    Kb_(Ru_/NA_)
{
    // determine the number is the specie which nucleats
    forAll(specieThermo_, specieI)
    {
        if (owner.slgThermo().carrier().species()[specieI] == nuclSpecies_)
        {
            nuclSpeciesIndex_ = specieI;
        }

        if (owner.slgThermo().carrier().species()[specieI] == inertSpecies_)
        {
            inertSpeciesIndex_ = specieI;
        }
    }
        
    molecularWeightNuclSpecies_ = specieThermo_[nuclSpeciesIndex_].W();
    molecularMass_ = molecularWeightNuclSpecies_/NA_;
    
    Info << "Nucleating species: " << owner.slgThermo().carrier().species()[nuclSpeciesIndex_] << endl;
    Info << "Inert species: " << owner.slgThermo().carrier().species()[inertSpeciesIndex_] << endl;
    Info << "Molecular weight of nucleating species: " << molecularWeightNuclSpecies_ << endl;
    Info << "Molecular mass of nucleating species: " << molecularMass_ << endl;
    Info << "Diameter of first section: " << d0_ << endl;
    Info << "Value of volumetric growth rate: " << volGrowthFactor_ << endl;
    Info << "Avogadro-Constant: " << NA_ << endl;
    Info << "Gas constant: " << Ru_ << endl;
    Info << "Boltzman constant: " << Kb_ << endl;
    Info << "Formula for surface tension: " << surfaceTensionModelCoeffM_ 
         << "T + " << surfaceTensionModelCoeffN_ << endl;
    Info << "Formula for liquid density: " << liquidDensityModelCoeffM_ 
         << "T + " << liquidDensityModelCoeffN_ << endl;
    Info << "Formula for saturation pressure model: 10^(" << saturationPressureModelCoeffA_ 
         << "T^2 + " << saturationPressureModelCoeffB_ 
         << "T + " << saturationPressureModelCoeffC_ << ")" << endl;
    Info << "Selected nucleation rate model: " << nucleationRateModel_ << endl;
    Info << "Selected condensation rate model: " << condensationRateModel_ << endl;
    
    setSectionProperties();
}


template <class CloudType>
Foam::SectionalHomogeneousNucleationCondensation<CloudType>::SectionalHomogeneousNucleationCondensation
(
    const SectionalHomogeneousNucleationCondensation<CloudType>& cm
)
:
    SynthesisModel<CloudType>(cm),
    
    specieThermo_(cm.specieThermo_),
    
    saturationIndex_(0),
    
    nucleationRateIndex_(1),
    
    condensationRateIndex_(2),
    
    nPhysicalAerosolProperties_(3),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
         
    nuclSpecies_(this->coeffDict().lookup("nucleatingSpecies")),

    inertSpecies_(this->coeffDict().lookup("inertSpecie")),
    
    surfaceTensionModelCoeffs_(this->coeffDict().subDict("surfaceTensionModelCoeffs")),

    surfaceTensionModelCoeffM_(readScalar(surfaceTensionModelCoeffs_.lookup("m"))),

    surfaceTensionModelCoeffN_(readScalar(surfaceTensionModelCoeffs_.lookup("n"))),

    liquidDensityModelCoeffs_(this->coeffDict().subDict("liquidDensityModelCoeffs")),

    liquidDensityModelCoeffM_(readScalar(liquidDensityModelCoeffs_.lookup("m"))),

    liquidDensityModelCoeffN_(readScalar(liquidDensityModelCoeffs_.lookup("n"))),

    saturationPressureModelCoeffs_(this->coeffDict().subDict("saturationPressureModelCoeffs")),

    saturationPressureModelCoeffA_(readScalar(saturationPressureModelCoeffs_.lookup("a"))),

    saturationPressureModelCoeffB_(readScalar(saturationPressureModelCoeffs_.lookup("b"))),

    saturationPressureModelCoeffC_(readScalar(saturationPressureModelCoeffs_.lookup("c"))),   
 
    nucleationRateModel_(this->coeffDict().lookup("nucleationRateModel")),
    
    condensationRateModel_(this->coeffDict().lookup("condensationRateModel")),
    
    d_(),
    
    a_(),
    
    v_(),
    
    nuclSpeciesIndex_(0),
    
    inertSpeciesIndex_(0),
    
    surfaceTension_(0),
    
    liquidDensity_(0),
    
    molecularWeightNuclSpecies_(0),
    
    molecularMass_(0),
    
    molecularVolume_(0),
    
    molecularSurface_(0),
    
    molecularWeightMixture_(0),
    
    saturationPressure_(0),
    
    vaporPressure_(0),
        
    mixtureDensity_(0),
    
    nucleationRate_(0),
    
    condensationRate_(0),
    
    meanFreePath_(0),
    
    d0_(readScalar(this->coeffDict().lookup("d0"))),
       
    volGrowthFactor_(readScalar(this->coeffDict().lookup("volGrowthFactor"))),
    
    criticalDiameter_(0),
    
    criticalVolume_(0),
        
    condensedMassFraction_(0),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    Kb_(Ru_/NA_)
{
    // determine the number is the specie which nucleats
    forAll(specieThermo_, specieI)
    {
        if (cm.slgThermo().carrier().species()[specieI] == nuclSpecies_)
        {
            nuclSpeciesIndex_ = specieI;
        }

        if (cm.slgThermo().carrier().species()[specieI] == inertSpecies_)
        {
            inertSpeciesIndex_ = specieI;
        }
    }

    molecularWeightNuclSpecies_ = specieThermo_[nuclSpeciesIndex_].W();
    molecularMass_ = molecularWeightNuclSpecies_/NA_;

    Info << "Nucleating species: " << cm.slgThermo().carrier().species()[nuclSpeciesIndex_] << endl;
    Info << "Inert species: " << cm.slgThermo().carrier().species()[inertSpeciesIndex_] << endl;
    Info << "Molecular weight of nucleating species: " << molecularWeightNuclSpecies_ << endl;
    Info << "Molecular mass of nucleating species: " << molecularMass_ << endl;
    Info << "Diameter of first section: " << d0_ << endl;
    Info << "Value of volumetric growth rate: " << volGrowthFactor_ << endl;
    Info << "Avogadro-Constant: " << NA_ << endl;
    Info << "Gas constant: " << Ru_ << endl;
    Info << "Boltzman constant: " << Kb_ << endl;
    Info << "Formula for surface tension: " << surfaceTensionModelCoeffM_ 
         << "T + " << surfaceTensionModelCoeffN_ << endl;
    Info << "Formula for liquid density: " << liquidDensityModelCoeffM_ 
         << "T + " << liquidDensityModelCoeffN_ << endl;
    Info << "Formula for saturation pressure model: 10^(" << saturationPressureModelCoeffA_ 
         << "T^2 + " << saturationPressureModelCoeffB_ 
         << "T + " << saturationPressureModelCoeffC_ << ")" << endl;
    Info << "Selected nucleation rate model: " << nucleationRateModel_ << endl;
    Info << "Selected condensation rate model: " << condensationRateModel_ << endl;
    
    setSectionProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SectionalHomogeneousNucleationCondensation<CloudType>::~SectionalHomogeneousNucleationCondensation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::synthesize()
{
    Info << "Nucleation and condensation..." << endl;    
    
    forAllIters(this->owner(), particle)
    {
        calcSurfaceTension(particle());
        calcDensityLiquid(particle());
        calcMolecularVolume();
        calcMolecularSurface();
        calcMolecularWeightMixture(particle());
        calcSaturationPressure(particle());
        calcVaporPressure(particle());
        calcSaturation(particle());
        calcCriticalDiameter(particle());
        calcCriticalVolume();
        calcMixtureDensity(particle());
        calcNucleationRate(particle());
        nucleation(particle());
        condensation(particle());
        latentHeat(particle());
        substractCondensedMass(particle());
    }
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::setSectionProperties()
{
    d_.setSize(nSections_, 0.0);
    a_.setSize(nSections_, 0.0);
    v_.setSize(nSections_, 0.0);

    Info << endl << "section" << token::TAB << "diameter" << token::TAB << "surface area" << token::TAB << "volume" << endl;
    
    v_[0] = M_PI/6.0*Foam::pow(d0_, 3.0);
    
    forAll (v_, section)
    {
        v_[section] = v_[0]*Foam::pow(volGrowthFactor_,section);
        d_[section] = Foam::pow((6.0*v_[section]/M_PI),1.0/3.0);
        a_[section] = M_PI*Foam::pow(d_[section], 2.0);
        
        Info << section << token::TAB << d_[section] << token::TAB << a_[section]  << token::TAB << v_[section] << endl;
    }
    Info << endl;
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcSurfaceTension(particleType& particle) // kg/s^2
{
    surfaceTension_ = surfaceTensionModelCoeffM_*particle.T() + surfaceTensionModelCoeffN_;
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcDensityLiquid(particleType& particle) // kg/m^3
{
    liquidDensity_ = liquidDensityModelCoeffM_*particle.T() + liquidDensityModelCoeffN_; 
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcMixtureDensity(particleType& particle) // kg/m^3
{
    mixtureDensity_ = (particle.pc()*molecularWeightMixture_)/(Ru_*particle.T()); // Warnatz: Combustion (2001), Equation (1.5), Page 4
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcMolecularVolume() // m^3
{
    molecularVolume_ = (molecularMass_/liquidDensity_); // Di Veroli, Rigopoulos: Physics of Fluids 23, 043305 (2011), Page 6
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcMolecularSurface() // m^2
{
    molecularSurface_ = (Foam::pow(4.0*M_PI,(1./3.))*Foam::pow(3.0*molecularVolume_,(2./3.)));
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcMolecularWeightMixture(particleType& particle) // g/mol
{
    molecularWeightMixture_ = 0.0;

    forAll(specieThermo_,specieI) 
    {
        molecularWeightMixture_ += particle.Y()[specieI]/specieThermo_[specieI].W();
    }

    molecularWeightMixture_ = 1.0/molecularWeightMixture_;
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcSaturationPressure(particleType& particle) // 1Pa = kg*m/s^2
{
    saturationPressure_ = Foam::pow(10., saturationPressureModelCoeffA_*Foam::sqr(particle.T()) + saturationPressureModelCoeffB_*particle.T() + saturationPressureModelCoeffC_);
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcVaporPressure(particleType& particle) 
{
    vaporPressure_ = (particle.Y()[nuclSpeciesIndex_]*molecularWeightMixture_/molecularWeightNuclSpecies_*particle.pc());     
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcSaturation(particleType& particle) // [-]            
{
    if (saturationPressure_ > VSMALL) 
        particle.physicalAerosolProperties()[saturationIndex_] = vaporPressure_/saturationPressure_;
    else 
        particle.physicalAerosolProperties()[saturationIndex_] =  0; 
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcCriticalDiameter(particleType& particle) // [m]
{
    if (particle.physicalAerosolProperties()[saturationIndex_] > SMALL)
        criticalDiameter_ = ((4*surfaceTension_*molecularVolume_)/(Kb_*particle.T()*Foam::log(particle.physicalAerosolProperties()[saturationIndex_]))); // Di Veroli, Rigopoulos: Physics of Fluids 23, 043305 (2011), Page 6
    else
        criticalDiameter_ = 0;
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcCriticalVolume() // [m^3]
{
    criticalVolume_ = 4./3.*M_PI*Foam::pow((criticalDiameter_/2.),3.);
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcMeanFreePath(particleType& particle) // [m]
{
    meanFreePath_ = (Kb_*particle.T()/(Foam::pow(2.0, 0.5)*M_PI*particle.pc()*Foam::pow(4.615e-10,2.0))); // Wikipedia
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcNucleationRate(particleType& particle) 
{
    nucleationRate_ = 0.0;
    
    if (particle.physicalAerosolProperties()[saturationIndex_] > SMALL) 
    {
        scalar partialDensity = particle.Y()[nuclSpeciesIndex_]*mixtureDensity_; // [kg/m^3]
        
        // get number of molecules/Volume (concentration of vapor)
        scalar Nv = partialDensity*NA_/molecularWeightNuclSpecies_; // [m^-3]
        
        const scalar baseFormulation = Foam::pow(Nv,2.)*molecularVolume_
                                     * Foam::pow((2.*surfaceTension_/(M_PI*molecularMass_)),0.5)
                                     * Foam::exp((-16.*M_PI*Foam::pow(surfaceTension_,3.)
                                     * Foam::pow(molecularMass_,2.))/(3.*Foam::pow(Kb_*particle.T(),3.)
                                     * Foam::pow(liquidDensity_,2.)
                                     * Foam::pow(Foam::log(particle.physicalAerosolProperties()[saturationIndex_]),2.)));
        
        if (nucleationRateModel_ == "gamory")
        {
            nucleationRate_ = baseFormulation;
        }
        else if (nucleationRateModel_ == "veroli")
        {            
            nucleationRate_ = baseFormulation*Foam::exp((molecularSurface_*surfaceTension_)/(Kb_*particle.T()));
        }
        else if (nucleationRateModel_ == "girshick")
        {
            nucleationRate_ = baseFormulation*Foam::exp((molecularSurface_*surfaceTension_)/(Kb_*particle.T()))*(1./particle.physicalAerosolProperties()[saturationIndex_]);
        }
        else if(nucleationRateModel_ == "none")
        {
            nucleationRate_ = 0.0;
        }
        else 
        {
            std::stringstream sstm;
            sstm << "     Foam::SectionalHomogeneousNucleationCondensation::calcNucleationRate() -- no nucleation model selected";
            sstm << "     Valid nucleation models are:" << endl;
            sstm << "                   girshick" << endl;
            sstm << "                   garmory" << endl;
            sstm << "                   veroli" << endl;
            sstm << "                   none" << endl;
            error myError(sstm.str());
            myError.exit(1);
        }
    }
    
    if (nucleationRate_ < 0)
    {
        nucleationRate_ = 0;
        
        Info << "Negativ nucleation!!!";
    }
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcCondensationRate(particleType& particle, scalar d)
{
    // diffusion coefficient [m2/s]
    const scalar Diff = 0.000025*(Foam::pow(particle.T(),1.75)/particle.pc());
    
    if (condensationRateModel_ == "freeMolecular")
    {
        condensationRate_ = ((2.0*molecularVolume_*saturationPressure_)/Foam::sqrt(2*M_PI*molecularMass_*Kb_*particle.T()))*(particle.physicalAerosolProperties()[saturationIndex_] - 1.0);
    }
    else if (condensationRateModel_ == "continuum")
    {
        condensationRate_ = ((4.0*Diff*molecularVolume_*saturationPressure_)/(Kb_*particle.T()*d))*(particle.physicalAerosolProperties()[saturationIndex_] - 1.0);
    }
    else if (condensationRateModel_ == "harmonicMean")
    {
        scalar condensationRatefreeMolecular = ((2.0*molecularVolume_*saturationPressure_)/Foam::sqrt(2*M_PI*molecularMass_*Kb_*particle.T()))*(particle.physicalAerosolProperties()[saturationIndex_] - 1.0);
        
        scalar condensationRateContinuum     = ((4.0*Diff*molecularVolume_*saturationPressure_)/(Kb_*particle.T()*d))*(particle.physicalAerosolProperties()[saturationIndex_] - 1.0);
        
        condensationRate_                    = condensationRatefreeMolecular*condensationRateContinuum/(condensationRatefreeMolecular + condensationRateContinuum);
    }
    else if (condensationRateModel_ == "FuchsSutugin")
    {
        calcMeanFreePath(particle);
        
        // Knudsen number
        const scalar Kn = 2*meanFreePath_/d;
    
        // blend function
        const scalar alphaKn = (1 + Kn)/(1 + 1.71*Kn + 1.333*Foam::pow(Kn,2.));

        condensationRate_ = alphaKn*((4.0*Diff*molecularVolume_*saturationPressure_)/(Kb_*particle.T()*d))*(particle.physicalAerosolProperties()[saturationIndex_] - 1.0);
    }    
    else if (condensationRateModel_ == "FuchsSutugin")
    {
        condensationRate_ = 0.0;
    }
    else 
    {
        std::stringstream sstm;
        sstm << "     Foam::SectionalHomogeneousNucleationCondensation<CloudType>::calcCondensationRate -- no condensation model selected";
        sstm << "     Valid condensation models are:" << endl;
        sstm << "                   freeMolecular" << endl;
        sstm << "                   continuum" << endl;
        sstm << "                   harmonicMean" << endl;
        sstm << "                   FuchsSutugin" << endl;
        sstm << "                   none" << endl;
        error myError(sstm.str());
        myError.exit(1);
    }
    
    if (condensationRate_ < 0)
        condensationRate_ = 0.0;
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::nucleation(particleType& particle)
{
    CloudType& cloud(this->owner());
    const scalar deltaT = cloud.mesh().time().deltaT().value(); // [s]
    
    condensedMassFraction_ = 0.0;
    
    if ((particle.physicalAerosolProperties()[nucleationRateIndex_] > SMALL) && (criticalDiameter_ > d_[0]))
    {
        label  section = -1;
        scalar inputSize = 0;
        label  lowerSection = 0;
        label  upperSection = 0;
        
        // find closest section above critical diameter
        while (inputSize == 0)
        {
            section++;
            
            if (d_[section] > criticalDiameter_)
            {
                inputSize = d_[section]; 
            }
        }
        
        upperSection = section;
        lowerSection = section - 1;
        
        if (lowerSection < 0)
        {
            lowerSection = 0;
        }
        
        const scalar numberNuclDrop = particle.physicalAerosolProperties()[nucleationRateIndex_]*deltaT; // [#/m^3] = [#/m^3 s] [s]
        
        const scalar numberNuclDropUpper = numberNuclDrop*(v_[upperSection] - criticalVolume_)/(v_[upperSection] - v_[lowerSection]);
        const scalar numberNuclDropLower = numberNuclDrop - numberNuclDropUpper;
                
        particle.psdProperties()[lowerSection] += numberNuclDropLower;
        particle.psdProperties()[upperSection] += numberNuclDropUpper;
        
        condensedMassFraction_ = (numberNuclDropLower*v_[lowerSection] + numberNuclDropUpper*v_[upperSection])*liquidDensity_/mixtureDensity_; // [-] = [#/m^3] [m^3] [kg/m^3] / [kg/m^3]
    }
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::condensation(particleType& particle)
{
    CloudType& cloud(this->owner());
    const scalar deltaT = cloud.mesh().time().deltaT().value(); // [s]
    
    const label nSec = particle.psdProperties().size();
    
    scalarList sectionChange(nSec, 0.0);
    
    forAll (particle.psdProperties(), section)
    {   
        const scalar pi = constant::mathematical::pi;
        const scalar d = d_[section];
        
        calcCondensationRate(particle, d);
        
        if (condensationRate_ > 0.0)
        {
            const scalar targetDiameter = d + condensationRate_*deltaT;
            scalar targetVolume   = 1.0/6.0*pi*Foam::pow((targetDiameter), 3.0);
            
            if (targetVolume > v_[nSec-1])
            {
                targetVolume = v_[nSec-1];
            }
            
            label targetSectionL = section;
            label targetSectionU = section + 1;
            
            forAll(particle.psdProperties(), searchSection)
            {
                if ((targetVolume >= v_[searchSection]) && (targetVolume < v_[searchSection+1]))
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
        }
        
        if (section == 10)
        {
            particle.physicalAerosolProperties()[condensationRateIndex_] = condensationRate_;
        }
    }
    
    forAll (particle.psdProperties(), section)
    {
        particle.psdProperties()[section] += sectionChange[section];
        
        condensedMassFraction_ += sectionChange[section]*v_[section]*liquidDensity_/mixtureDensity_;
    }
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::latentHeat(particleType& particle)
{ 
    if (condensedMassFraction_ > VSMALL)
    {
        // heat of vaporization: 79.2 kJ/mol at 340 deg C 
        // Source: Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 95th Edition. CRC Press LLC, Boca Raton: FL 2014-2015, p. 6-135
        const scalar latentHeat = 79200000.0/molecularWeightNuclSpecies_; // [J/kg] = [J/kmol] / [kg/kmol]

        const scalar deltaH = condensedMassFraction_*latentHeat; // [J] = [J/kg] [kg]
        
        particle.hA() += deltaH/particle.wt(); // [J/kg] = [J] / [kg]
    }
}


template <class CloudType>
void Foam::SectionalHomogeneousNucleationCondensation<CloudType>::substractCondensedMass(particleType& particle)
{
    if (condensedMassFraction_ > VSMALL)
    {
        particle.Y()[nuclSpeciesIndex_]  = max(particle.Y()[nuclSpeciesIndex_] - condensedMassFraction_, 0.0); 
        particle.Y()[inertSpeciesIndex_] = min(particle.Y()[inertSpeciesIndex_] + condensedMassFraction_, 1.0); ;
        
        // reset condensed mass fraction after substraction
        condensedMassFraction_ = 0.0;
    }
}


template<class CloudType>
const Foam::wordList Foam::SectionalHomogeneousNucleationCondensation<CloudType>::psdPropertyNames() const
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
const Foam::wordList Foam::SectionalHomogeneousNucleationCondensation<CloudType>::physicalAerosolPropertyNames() const
{
    wordList physicalAerosolPropertyNames;
    physicalAerosolPropertyNames.setSize(nPhysicalAerosolProperties_);
    
    physicalAerosolPropertyNames[0] = "saturation";
    physicalAerosolPropertyNames[1] = "nucleationRate";
    physicalAerosolPropertyNames[2] = "condensationRate";
        
    return physicalAerosolPropertyNames;
}


// ************************************************************************* //

