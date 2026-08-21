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

#include "PremixedMixingPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::PremixedMixingPopeParticle<ParticleType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::setCellValues(cloud, td, dt, cellI);
    
    tetIndices tetIs = this->currentTetIndices();
    
    if (cloud.mixing().XiRNames().size())
    {
        const mmcVarSet& setOfXi(cloud.mixing().XiR());
        
        const labelHashTable& indexInXi  = setOfXi.rVarInXi();
        const labelHashTable& indexInXiR = setOfXi.rVarInXiR();
        
        forAllConstIter(wordList, cloud.mixing().XiRNames(), iter)
        {
        
            //- This part requires review 
            //- (For now it works fine just for interpolated) 
            
            if (XiR().size() > 1) 
            {
                FatalError << "Reference variable relaxation is implemented only for a SINGLE reference variable." 
                           << exit(FatalError);
            }
              
            if (setOfXi.Vars(indexInXi[*iter]).refType()=="interpolated")
            {
                progVarInt() = td.XiRInterp()[indexInXiR[*iter]].interpolate(this->coordinates(),tetIs);
                if (progVarInt() > 1.0) progVarInt() = 1.0;   // correction for odd interpolation results (e.g. on periodic b.c.)
                XiR()[indexInXiR[*iter]] = XiR()[indexInXiR[*iter]] + XiRrlxFactor() * (progVarInt() - XiR()[indexInXiR[*iter]]);   
            }             
            else if (setOfXi.Vars(indexInXi[*iter]).refType()=="evolved")
            {
                // Do nothing since it is evolved by MMC model
            }
            else
                XiR()[indexInXiR[*iter]] = 0.0;//this->position();
        }
    }
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::PremixedMixingPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::calc(cloud, td, dt, cellI);
}

// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::PremixedMixingPopeParticle<ParticleType>::mixProperties
(
    PremixedMixingPopeParticle<ParticleType>& p, 
    PremixedMixingPopeParticle<ParticleType>& q, 
    scalar mixExtent
)
{       
    ParticleType::mixProperties(p, q, mixExtent);
}


template<class ParticleType>
void Foam::PremixedMixingPopeParticle<ParticleType>::mixProperties
(
    PremixedMixingPopeParticle<ParticleType>& p, 
    PremixedMixingPopeParticle<ParticleType>& q, 
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{       
    ParticleType::mixProperties(p, q, mixExtent,mixExtentSoot);
}


template<class ParticleType>
void Foam::PremixedMixingPopeParticle<ParticleType>::mixProperties
(
    PremixedMixingPopeParticle<ParticleType>& p, 
    PremixedMixingPopeParticle<ParticleType>& q, 
    scalar mixExtent,
    scalarList ScaledExtent
)
{       
    ParticleType::mixProperties(p, q, mixExtent,ScaledExtent);
}


template<class ParticleType>
template<class CloudType>
void Foam::PremixedMixingPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);
    
    if (c.mixing().numXiR())
    {
        indexInXiR_ = c.mixing().XiR().rVarInXiR();
        XiRNames_ = c.mixing().XiRNames();
    }
}

// * * * * * * * * * *  Statistical Data Functions * * * * * * * * * * * * * * *
template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::PremixedMixingPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::PremixedMixingPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::PremixedMixingPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();



    // Add ethalpy
    this->nameVariableLookUpTable().addNamedVariable("dx",dx_);

    this->nameVariableLookUpTable().addNamedVariable("mixTime",mixTime_);
    
    this->nameVariableLookUpTable().addNamedVariable("mixExt",mixExt_);
    
    this->nameVariableLookUpTable().addNamedVariable("progVarInt",progVarInt_);

    this->nameVariableLookUpTable().addNamedVariable("XiRrlxFactor",XiRrlxFactor_);
    
    this->nameVariableLookUpTable().addNamedVariable("rlxStep",rlxStep_);

    // Note: XiRNames is initialized as a static variablebefore 
    // XiR_ is set. Therefore we need to check if they are set to avoid
    // accessing elements out of bounds

    // Add reference variables
    if (XiRNames_.size() == XiR_.size())
    {
        for (const word& name : XiRNames_)
        {
            this->nameVariableLookUpTable().addNamedVariable
            (
                name,
                XiR(name)
            );

            word dXirName = "d"+name;
            this->nameVariableLookUpTable().addNamedVariable
            (
                dXirName,
                dXiR(name)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::PremixedMixingPopeParticle<ParticleType>::PremixedMixingPopeParticle
(
    const PremixedMixingPopeParticle<ParticleType>& p
)
:
    ParticleType(p),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    mixTime_(p.mixTime_),
    mixExt_(p.mixExt_),
    progVarInt_(p.progVarInt_),
    XiRrlxFactor_(p.XiRrlxFactor_),
    rlxStep_(p.rlxStep_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::PremixedMixingPopeParticle<ParticleType>::PremixedMixingPopeParticle
(
    const PremixedMixingPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
:
    ParticleType(p, mesh),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    mixTime_(p.mixTime_),
    mixExt_(p.mixExt_),
    progVarInt_(p.progVarInt_),
    XiRrlxFactor_(p.XiRrlxFactor_),
    rlxStep_(p.rlxStep_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "PremixedMixingPopeParticleIO.C"


// ************************************************************************* //
