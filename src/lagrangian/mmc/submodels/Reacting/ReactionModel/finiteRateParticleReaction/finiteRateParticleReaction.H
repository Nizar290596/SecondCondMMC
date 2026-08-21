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
    Foam::finiteRateParticleReaction

Description

SourceFiles
    finiteRateParticleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef finiteRateParticleReaction_H
#define finiteRateParticleReaction_H

#include "BalanceReactModel.H"

#include "particleTDACChemistryModel.H"
#include "ode.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class finteRateParticleReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class finiteRateParticleReaction
:
    public ode
           <
                particleTDACChemistryModel
                <
                    psiReactionThermo,
                    gasHThermoPhysics
                >
           >,
           
    public BalanceReactModel<CloudType>
{ 
    private:
        scalar zLower_;
        scalar zUpper_;

    public:

        //- Runtime type information
        TypeName("finiteRateParticleChemistry");


    // Constructors

        //- Construct from dictionary
        finiteRateParticleReaction(const dictionary& dict, CloudType& owner);


    // Destructor
    virtual ~finiteRateParticleReaction();

    
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
            scalarField& pY
        ) override;

        virtual void ISATCntIncrmnt();

        virtual void resetCPUTimeAnalysis();

        virtual void writeCPUTimeAnalysis();

        virtual label nActiveSp();

        virtual void Sreact();

        virtual scalar getCPUtime();


};


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "finiteRateParticleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
