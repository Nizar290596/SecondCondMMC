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

#ifndef SectionalDropletSynthesis_H
#define SectionalDropletSynthesis_H

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
                    Class SectionalDropletSynthesis Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SectionalDropletSynthesis
:
    public SynthesisModel<CloudType>
{
    // Private data
        //- ThermoType  instantiated for
        typedef gasHThermoPhysics ThermoType;
    
        const PtrList<ThermoType>& speciesThermo_;
        
        typedef typename CloudType::particleType particleType; 
        
        // physical aerosol properties
        const label saturation_;
        const label nucleationRate_;
        const label surfaceGrowthRate_;
        
        const label nPhysicalAerosolProperties_;
        
        const label nSections_;
        
        const scalar d0_;        
        const scalar growth_;
        
        // diameter, surface area and volume of the sections
        scalarList d_;
        scalarList a_;
        scalarList v_;
        
        scalarList X_;
        
        // species names
        const word condensingSpeciesName_;
        
        // indexes of important species
        label condensingSpeciesIndex_;
        
        // density of the solid species
        scalar liquidDensity_;
        
        // coagulation model
        const word coagulationModel_;
        
        // cool rate
        scalar coolrate_;
        
        // Avogadro constant
        scalar NA_;

        // universal gas constant
        scalar Ru_;
        
        // Boltzmann constant
        scalar kb_;
        
    public:

        //- Runtime type information
        TypeName("sectionalDropletSynthesis");

    // Constructors

        //- Construct from dictionary
        SectionalDropletSynthesis(const dictionary& dict, CloudType& owner);
        
        //- Construct copy
        SectionalDropletSynthesis(const SectionalDropletSynthesis<CloudType>& cm);

        //- Destructor
        virtual ~SectionalDropletSynthesis();

    // Member functions
        
        virtual void synthesize();
        
        virtual const wordList psdPropertyNames() const;
        
        virtual const wordList physicalAerosolPropertyNames() const;
        
        void setSectionProperties();
        
        void setSizeSplittingOperator();
        
        void particleInception(particleType& particle);
        
        void particleSurfaceGrowth(particleType& particle);
        
        void particleCoagulation(particleType& particle);
        
        void gasPhaseCooling(particleType& particle);
      
        scalar surfaceTension(scalar T);
        
        scalar dynamicViscosity(scalar T);
        
        scalar saturationPressure(scalar T);
        
        scalar calcMolecularWeightMixture(scalarField& Y);
        
        scalar meanFreePath(scalar T, scalar p, scalar mu);
        
        scalar criticalNode(scalar sigma, scalar T, scalar S);
        
        scalar criticalVolume(scalar sigma, scalar T, scalar S);
        
        scalar nucleationRate(scalar sigma, scalar T, scalar S, scalar ns);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SectionalDropletSynthesis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
