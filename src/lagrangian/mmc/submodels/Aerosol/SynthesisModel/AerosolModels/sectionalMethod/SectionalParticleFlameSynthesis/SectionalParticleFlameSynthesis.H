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

  Class
  Foam::NoReaction

  Description
  Dummy class for 'none' option - will raise an error if any functions are
  called that require return values.

  SourceFiles
  NoReaction.C

  \*---------------------------------------------------------------------------*/

#ifndef SectionalParticleFlameSynthesis_H
#define SectionalParticleFlameSynthesis_H

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
                    Class SectionalParticleFlameSynthesis Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SectionalParticleFlameSynthesis
:
    public SynthesisModel<CloudType>
{
    // Private data
        //- ThermoType  instantiated for
        typedef gasHThermoPhysics ThermoType;
    
        const PtrList<ThermoType>& speciesThermo_;
        
        typedef typename CloudType::particleType particleType;
        
        // physical aerosol properties
        const label inceptionRateIndex_;
        const label surfaceGrowthRateIndex_;
        const label primaryParticleNumberIndex_;
        const label particleNumberIndex_;
        const label geoMeanParticleVolumeIndex_;
        const label geoStdDevParticleVolumeIndex_;
        const label elsGasIndex_;
        const label elsParticleIndex_;
        const label elsIndex_;
        
        const label nPhysicalAerosolProperties_;
        
        const label nSections_;
        
        // Switches for particle inception, surface growth and coagulation
        const Switch inception_;
        const Switch surfaceGrowth_;
        const Switch coagulation_;
        
        const scalar d00_;
        const scalar growth_;
        const scalar kf_;
        const scalar Df0_;
        const scalar coalescenceDiameter_;
        const scalar dDfLimit_;
        const scalar DfLimit_;
        const scalar growthLimit_;
	const scalar elsGasNorm_;
        const scalar elsParticleNorm_;
        const scalar zLower_;
        const scalar zUpper_;
        
        // diameter, surface area and volume of the sections
        scalarList d_;
        scalarList a_;
        scalarList v_;
        
        scalarList d0_;
        scalarList dc_;
        scalarList N_;
        
        scalarList Df_;
        scalarList Rg_;
        
        scalarList X_;
        scalarList beta_;
        
        // species names
        const dictionary precursorSpeciesDict_;
        const word solidSpeciesName_;
                
        // indexes of important species
        label precursorSpeciesIndex_;
        label depositionSpeciesIndex_;
        label solidSpeciesIndex_;        
                
        // density of the solid species
        const scalar solidDensity_;
        
        // surface reaction model
        const scalar surfaceReactionEfficiency_;
        const dictionary depositingSpeciesDict_;        
        
        // coagulation model
        const word coagulationModel_;
        
        const Switch turbCoagulation_;

        const label coagulationSubSteps_;
        
        // fields for turbulent coagulation
        const volScalarField& mu_;
        const volScalarField& nut_; 
        const volScalarField& epsilon_;
        
        // Avogadro constant
        const scalar NA_;

        // universal gas constant
        const scalar Ru_;
        
        // Boltzmann constant
        const scalar kB_;
        
    public:

        //- Runtime type information
        TypeName("sectionalParticleFlameSynthesis");

    // Constructors

        //- Construct from dictionary
        SectionalParticleFlameSynthesis(const dictionary& dict, CloudType& owner);
        
        //- Construct copy
        SectionalParticleFlameSynthesis(const SectionalParticleFlameSynthesis<CloudType>& cm);

        //- Destructor
        virtual ~SectionalParticleFlameSynthesis();

    // Member functions
        
        virtual void synthesize();
        
        virtual const wordList psdPropertyNames() const;
        
        virtual const wordList physicalAerosolPropertyNames() const;
        
        void setSectionProperties();
        
        void setSizeSplittingOperator();
        
        void particleInception(particleType& particle);
        
        void particleSurfaceGrowth(particleType& particle);
        
        void particleCoagulation(particleType& particle);
        
        void calcPhysicalAerosolProperties(particleType& particle);
        
        scalar calcMolecularWeightMixture(scalarField& Y);

        scalarList X(scalarField& Y);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SectionalParticleFlameSynthesis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
