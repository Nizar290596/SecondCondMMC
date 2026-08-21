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

#include "sootExtern.H"
#include "sootParticleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::sootParticleReaction<CloudType>::sootParticleReaction
(
    const dictionary& dict,
    CloudType& owner
)
:
    ode
    <
        particleChemistryModel
        <
            psiReactionThermo,
            gasHThermoPhysics
        >
    >(owner.thermo()),
    
    ReactionModel<CloudType>(dict,owner,typeName),

    zLower_(this->coeffDict().lookupOrDefault("zLower", 0.0)),
    
    zUpper_(this->coeffDict().lookupOrDefault("zUpper", 1.0))
{
    startup_
    (
    );
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::sootParticleReaction<CloudType>::~sootParticleReaction()
{}

// ************************************************************************* //

template<class CloudType>
void Foam::sootParticleReaction<CloudType>::calculate
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
    if (pz < zLower_ || pz > zUpper_)
    {
        return;
    }

    label nspec = Y.size();
    
    scalar spec0[nspec];
    
    forAll(Y,np)
    {
        spec0[np] = Y[np];
    }
    
    scalar Temp = pT;
    
    scalar pascalTokPa = 0.001;
    scalar Press = pc * pascalTokPa;
    
    particledatainput_
    (
        spec0,
        &Temp,
        &Press
    );
    
    scalar delt = dt;
    eachparticle_
    (
        &delt,
        spec0,
        &Temp,
        &Press
    );

    scalar spect[nspec];
    
    particledataoutput_
    (
        spect
    );
    
    forAll(Y,np)
    {
        Y[np] = spect[np];
    } 

}


// ************************************************************************* //

