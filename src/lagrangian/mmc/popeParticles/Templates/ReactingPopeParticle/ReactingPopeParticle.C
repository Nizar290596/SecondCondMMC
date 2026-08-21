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

#include "ReactingPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::ReactingPopeParticle<ParticleType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::setCellValues(cloud, td, dt, cellI);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::ReactingPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::calc(cloud, td, dt, cellI);

    // Non-subset particles skip chemistry integration.
    // Transport (position SDE) and XiR evolution are handled by ParticleType::calc() above.
    if (cloud.secondCondMixingEnabled() && this->secondCondFlag() == 0)
        return;

//    Info << "calc reaction" << endl;
    
    // Reaction models uses mixture fraction. For now I pass the conditioning 
    // variable  used for density coupling but something more generic is required.

    // Conditioning (state) variable
    const word cVarName = cloud.coupling().cVarName();

    if(!cloud.balanceReactionLoad())
    {
        scalar t0 = cloud.db().time().value();

        cloud.reaction().calculate
        (
            t0,dt,this->hA(),this->pc(),this->T(),this->XiC(cVarName),this->Y()
        );
    }


    this->T() = cloud.composition().particleMixture
    (
        this->Y()
    ).THa(this->hA(), this->pc(), this->T()
    );
//    Info << "T after update from raction ="<< this->T() << endl;
    this-> NumActSp() = cloud.reaction().nActiveSp();
}


template<class ParticleType>
void Foam::ReactingPopeParticle<ParticleType>::mixProperties
(
    ReactingPopeParticle<ParticleType>& p, 
    ReactingPopeParticle<ParticleType>& q,
    scalar mixExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent);
}


template<class ParticleType>
void Foam::ReactingPopeParticle<ParticleType>::mixProperties
(
    ReactingPopeParticle<ParticleType>& p, 
    ReactingPopeParticle<ParticleType>& q,
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{
    ParticleType::mixProperties(p, q, mixExtent, mixExtentSoot);
}


template<class ParticleType>
void Foam::ReactingPopeParticle<ParticleType>::mixProperties
(
    ReactingPopeParticle<ParticleType>& p, 
    ReactingPopeParticle<ParticleType>& q, 
    scalar mixExtent,
    scalarList ScaledExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent,ScaledExtent);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReactingPopeParticle<ParticleType>::ReactingPopeParticle
(
    const ReactingPopeParticle<ParticleType>& p
)
    :
    ParticleType(p)
{}


template<class ParticleType>
Foam::ReactingPopeParticle<ParticleType>::ReactingPopeParticle
(
    const ReactingPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
    :
    ParticleType(p, mesh)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingPopeParticleIO.C"

// ************************************************************************* //

