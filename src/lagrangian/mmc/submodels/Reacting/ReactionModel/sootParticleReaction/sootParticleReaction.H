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
  Foam::sootParticleReaction

  Description

  SourceFiles
  sootParticleReaction.C

  \*---------------------------------------------------------------------------*/

#ifndef sootParticleReaction_H
#define sootParticleReaction_H

#include "ReactionModel.H"

#include "particleChemistryModel.H"
#include "ode.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class sootParticleReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class sootParticleReaction
:
    public ode
           <
                particleChemistryModel
                <
                    psiReactionThermo,
                    gasHThermoPhysics
                >
           >,
           
    public ReactionModel<CloudType>
{
    
    // Private data
    
        //- Upper and lower mixture fraction limits for which chemical reaction rates are calculated
        const scalar zLower_;
        const scalar zUpper_;
        
    public:

        //- Runtime type information
        TypeName("sootParticleChemistry");


    // Constructors

        //- Construct from dictionary
        sootParticleReaction(const dictionary& dict, CloudType& owner);


    // Destructor
    virtual ~sootParticleReaction();

    
    // Member functions
    
        //- Calculate chemistry
        virtual void calculate
        (
            const scalar t0, 
            const scalar dt,
            const scalar ha, 
            const scalar pc, 
            const scalar pT,
            const scalar pz,
            scalarField& Y
        ) override;
        
};


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sootParticleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
