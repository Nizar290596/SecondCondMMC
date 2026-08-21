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

#include "ReactingPopeCloud.H"
#include "ReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::setModels()
{
    reactionModel_.reset
    (
        ReactionModel<ReactingPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

}


template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::cloudReset(ReactingPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingPopeCloud<CloudType>::ReactingPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        U,
        DEff,
        rho,
        gradRho,
        Xi,
        false,
        true
    ),

    reactingPopeCloud(),

    cloudCopyPtr_(nullptr),

    balanceReactionLoad_(false),

    reactionModel_(nullptr)
{
    Info << "Creating reacting Pope Particle Cloud." << nl << endl;

    setModels();
    
    Info << nl << "Reaction model constructed." << endl;

    balanceReactionLoad_ = 
        this->subModelProperties().lookupOrDefault("balanceReactionLoad",false);

    Info << nl << "balanceReactionLoad is " << balanceReactionLoad() << endl;

    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Reacting Pope particle cloud data from file." << endl;
        
            particleType::reactingParticleIOType::readFields(*this);
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of reacting Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
                
                this->reaction().ignite();
            }
        }
    }
}


template<class CloudType>
Foam::ReactingPopeCloud<CloudType>::ReactingPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        Xi,
        false,      // Only top level cloud calls initAtCnstr
        readFields
    ),

    reactingPopeCloud(),

    cloudCopyPtr_(nullptr),

    balanceReactionLoad_(false),

    reactionModel_(nullptr)
{
    Info << "Creating reacting Pope Particle Cloud." << nl << endl;

    setModels();
    
    Info << nl << "Reaction model constructed." << endl;

    balanceReactionLoad_ = 
        this->subModelProperties().lookupOrDefault("balanceReactionLoad",false);

    Info << nl << "balanceReactionLoad is " << balanceReactionLoad() << endl;

    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Reacting Pope particle cloud data from file." << endl;
        
            particleType::reactingParticleIOType::readFields(*this);
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of reacting Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
                
                this->reaction().ignite();
            }
        }
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingPopeCloud<CloudType>::~ReactingPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle, mass, wt, patchI, patchFace,iniRls);
}


template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::setEulerianStatistics()
{}


template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();
}


template<class CloudType>
void Foam::ReactingPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::reactingParticleIOType::writeFields(*this);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::ReactingPopeCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename particleType::trackingData& td
)
{
    // when LoadBalance is on. solve does not compute chemistry.
    // so the time analysis are all 0
    
    if(!balanceReactionLoad())
        reactionModel_->resetCPUTimeAnalysis();

    CloudType::solve(cloud,td);
    
    if(!balanceReactionLoad())
    {
        reactionModel_->writeCPUTimeAnalysis();
        reactionModel_->ISATCntIncrmnt();
    }
}

// ************************************************************************* //
