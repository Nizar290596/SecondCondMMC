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

Class
    Foam::SectionalHomogeneousNucleationCondensation

Description

SourceFiles
    SectionalHomogeneousNucleationCondensation.C

\*---------------------------------------------------------------------------*/

#ifndef SectionalHomogeneousNucleationCondensation_H
#define SectionalHomogeneousNucleationCondensation_H

#include "SynthesisModel.H"

#include "List.H"
#include "scalar.H"
#include "scalarField.H"
#include "physicoChemicalConstants.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class SectionalHomogeneousNucleationCondensation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SectionalHomogeneousNucleationCondensation
:
    public SynthesisModel<CloudType>
{
    // Private data
  
        //- ThermoType instantiated for
        typedef gasHThermoPhysics ThermoType;
    
        const PtrList<ThermoType>& specieThermo_;
        
        typedef typename CloudType::particleType particleType; 
        
                // physical aerosol properties
        const label saturationIndex_;
        const label nucleationRateIndex_;
        const label condensationRateIndex_;
        
        const label nPhysicalAerosolProperties_;
        
        const label nSections_;
        
        // model selection
        const word nuclSpecies_;
        const word inertSpecies_;

        // surface tension model
        const dictionary surfaceTensionModelCoeffs_;
        const scalar surfaceTensionModelCoeffM_;
        const scalar surfaceTensionModelCoeffN_;

        // surface tension model
        const dictionary liquidDensityModelCoeffs_;
        const scalar liquidDensityModelCoeffM_;
        const scalar liquidDensityModelCoeffN_;

        // surface tension model
        const dictionary saturationPressureModelCoeffs_;
        const scalar saturationPressureModelCoeffA_;
        const scalar saturationPressureModelCoeffB_;
        const scalar saturationPressureModelCoeffC_;
        
        //- nucleation rate model
        const word nucleationRateModel_;
        
        //- condensation rate model
        const word condensationRateModel_;
        
        scalarList d_;
        scalarList a_;
        scalarList v_;
        
        // index of nucleating species
        label nuclSpeciesIndex_;
        
        // index of inert species
        label inertSpeciesIndex_;
        
        // surface tension
        scalar surfaceTension_;
        
        // density of liquid mass
        scalar liquidDensity_;
        
        // mulecular weight of the nucleation specie
        scalar molecularWeightNuclSpecies_;
        
        // mulecular mass of nucleation species
        scalar molecularMass_;
        
        // mulecular volume of the nucleation species
        scalar molecularVolume_;
        
        // mulecular surface of the nucleation species
        scalar molecularSurface_;
        
        // molecular mass of mixture
        scalar molecularWeightMixture_;
        
        // saturation pressure
        scalar saturationPressure_;

        // vapor pressure
        scalar vaporPressure_;
        
        // density of the gas mixture
        scalar mixtureDensity_;
        
        // nucleation rate
        scalar nucleationRate_;
        
        // condensation rate
        scalar condensationRate_;
        
        // mean free path of nucleating species molecules in inert species
        scalar meanFreePath_;
        
        // diameter of the primary particles
        scalar d0_;
        
        // volumetric growth factor
        scalar volGrowthFactor_;
        
        // critical particle diameter for nucleation
        scalar criticalDiameter_;
        
        // critical particle volume for nucleation
        scalar criticalVolume_;
        
        // accumulated nucleated mass which need to be subtracted from the gas phase
        scalar condensedMassFraction_;
        
        // Avogadro constant
        const scalar NA_;

        // universal gas constant
        const scalar Ru_;
        
        // Boltzmann constant
        const scalar Kb_;
                      
    public:

        //- Runtime type information
        TypeName("sectionalHomogeneousNucleationCondensation");

    // Constructors

        //- Construct from dictionary
        SectionalHomogeneousNucleationCondensation(const dictionary& dict, CloudType& owner);
        
        //- Construct copy
        SectionalHomogeneousNucleationCondensation(const SectionalHomogeneousNucleationCondensation<CloudType>& cm);

        //- Destructor
        virtual ~SectionalHomogeneousNucleationCondensation();

    // Member functions
        
        virtual void synthesize();
        
        virtual const wordList psdPropertyNames() const;
        
        virtual const wordList physicalAerosolPropertyNames() const;
        
        void setSectionProperties();
        
        void calcSurfaceTension(particleType& particle);
        void calcDensityLiquid(particleType& particle);
        void calcMolecularVolume();
        void calcMolecularSurface();
        void calcMolecularWeightMixture(particleType& particle);
        void calcSaturationPressure(particleType& particle);
        void calcVaporPressure(particleType& particle);
        void calcSaturation(particleType& particle);
        void calcCriticalDiameter(particleType& particle);
        void calcCriticalVolume();
        void calcMeanFreePath(particleType& particle);
        void calcMixtureDensity(particleType& particle);
        void calcNucleationRate(particleType& particle);
        void nucleation(particleType& particle);
        void condensation(particleType& particle);
        void latentHeat(particleType& particle);
        void substractCondensedMass(particleType& particle);
        
        void calcCondensationRate
        (
            particleType& particle,
            scalar d
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SectionalHomogeneousNucleationCondensation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
