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

#include "finiteRateParticleReaction.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::finiteRateParticleReaction<CloudType>::finiteRateParticleReaction
(
    const dictionary& dict,
    CloudType& owner
)
:
    ode
    <
        particleTDACChemistryModel
        <
            psiReactionThermo,
            gasHThermoPhysics
        >
    >(dynamic_cast<psiReactionThermo&>(owner.thermo())),
    
    BalanceReactModel<CloudType>(dict,owner,typeName),
    
    zLower_(this->coeffDict().getOrDefault("zLower", 0.0)),
    
    zUpper_(this->coeffDict().getOrDefault("zUpper", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::finiteRateParticleReaction<CloudType>::~finiteRateParticleReaction()
{}

// ************************************************************************* //

template<class CloudType>
void Foam::finiteRateParticleReaction<CloudType>::calculate
(
    const scalar t0, 
    const scalar dt, 
    const scalar ha,
    const scalar pc,
    const scalar pT,
    const scalar pz,
    scalarField& Y
)
{
    if (pz >= zLower_ && pz <= zUpper_)
    {
        particleCalculate(t0,dt,ha,pc,pT,Y);
    }
}

template<class CloudType>
void Foam::finiteRateParticleReaction<CloudType>::ISATCntIncrmnt()
{
    Info<<"\tUsing finiteRate Chemistry"<<endl;
    this->incrsTimeStepCounter(); 
    //timestep used by ISAT to clean tree etc.
    
}

template<class CloudType>
void Foam::finiteRateParticleReaction<CloudType>::resetCPUTimeAnalysis()
{
    this->setZeroCPUTimeData();
}

template<class CloudType>
void Foam::finiteRateParticleReaction<CloudType>::writeCPUTimeAnalysis()
{
    this->writeCPUTimeLog();
}

template<class CloudType>
Foam::label Foam::finiteRateParticleReaction<CloudType>::nActiveSp()
{
    label numActiveSpecies_(0);

    if(mechRed()->active())
    numActiveSpecies_ = mechRed()->NsSimp();

    return numActiveSpecies_;
}

template<class CloudType>
void Foam::finiteRateParticleReaction<CloudType>::Sreact()
{
    this->SreactBalance();
}

template<class CloudType>
Foam::scalar Foam::finiteRateParticleReaction<CloudType>::getCPUtime()
{
    return getCPUTimeReaction();
}


// ************************************************************************* //

