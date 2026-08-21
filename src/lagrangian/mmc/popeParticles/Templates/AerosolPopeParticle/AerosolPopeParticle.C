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

#include "AerosolPopeParticle.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::AerosolPopeParticle<ParticleType>::setCellValues
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
void Foam::AerosolPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::calc(cloud, td, dt, cellI);
}


template<class ParticleType>
void Foam::AerosolPopeParticle<ParticleType>::mixProperties
(
    AerosolPopeParticle<ParticleType>& p, 
    AerosolPopeParticle<ParticleType>& q, 
    scalar mixExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent);
    
    // mix particle PSD properties

    forAll(p.psdProperties(),n)
    { 
        //- Weighted pair mean
        scalar XiCAv = (p.wt() * p.psdProperties()[n] + q.wt() * q.psdProperties()[n])/(p.wt() + q.wt());

        //- mix them
        p.psdProperties()[n] = p.psdProperties()[n] + mixExtent * (XiCAv - p.psdProperties()[n]);
        q.psdProperties()[n] = q.psdProperties()[n] + mixExtent * (XiCAv - q.psdProperties()[n]);
    }
}


template<class ParticleType>
void Foam::AerosolPopeParticle<ParticleType>::mixProperties
(
    AerosolPopeParticle<ParticleType>& p, 
    AerosolPopeParticle<ParticleType>& q, 
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{
    ParticleType::mixProperties(p, q, mixExtent, mixExtentSoot);
    
    // mix particle PSD properties

    forAll(p.psdProperties(),n)
    { 
        //- Weighted pair mean
        scalar XiCAv = (p.wt() * p.psdProperties()[n] + q.wt() * q.psdProperties()[n])/(p.wt() + q.wt());

        //- mix them
        p.psdProperties()[n] = p.psdProperties()[n] + mixExtent * (XiCAv - p.psdProperties()[n]);
        q.psdProperties()[n] = q.psdProperties()[n] + mixExtent * (XiCAv - q.psdProperties()[n]);
    }
}


template<class ParticleType>
void Foam::AerosolPopeParticle<ParticleType>::mixProperties
(
    AerosolPopeParticle<ParticleType>& p, 
    AerosolPopeParticle<ParticleType>& q, 
    scalar mixExtent,
    scalarList ScaledExtent
)
{
    ParticleType::mixProperties(p, q, mixExtent, ScaledExtent);
    
    // mix particle PSD properties

    forAll(p.psdProperties(),n)
    { 
        //- Weighted pair mean
        scalar XiCAv = (p.wt() * p.psdProperties()[n] + q.wt() * q.psdProperties()[n])/(p.wt() + q.wt());

        //- mix them
        p.psdProperties()[n] = p.psdProperties()[n] + mixExtent * (XiCAv - p.psdProperties()[n]);
        q.psdProperties()[n] = q.psdProperties()[n] + mixExtent * (XiCAv - q.psdProperties()[n]);
    }
}


template<class ParticleType>
template<class CloudType>
void Foam::AerosolPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);
    
    physicalAerosolPropertyNames_ = c.synthesis().physicalAerosolPropertyNames();
    psdPropertyNames_ = c.synthesis().psdPropertyNames();
}



template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::AerosolPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::AerosolPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::AerosolPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();

    
    // Check that the static properties have already been initialized
    if (physicalAerosolPropertyNames_.size() == physicalAerosolProperties_.size())
    {
        forAll(physicalAerosolPropertyNames_,i)
        {
            this->nameVariableLookUpTable().addNamedVariable
            (
                physicalAerosolPropertyNames_[i],
                physicalAerosolProperties_[i]
            );
        }

        forAll(psdPropertyNames_,i)
        {
            this->nameVariableLookUpTable().addNamedVariable
            (
                psdPropertyNames_[i],
                psdProperties_[i]
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::AerosolPopeParticle<ParticleType>::AerosolPopeParticle
(
    const AerosolPopeParticle<ParticleType>& p
)
    :
    ParticleType(p),
    physicalAerosolProperties_(p.physicalAerosolProperties_),
    psdProperties_(p.psdProperties_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::AerosolPopeParticle<ParticleType>::AerosolPopeParticle
(
    const AerosolPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
    :
    ParticleType(p, mesh),
    physicalAerosolProperties_(p.physicalAerosolProperties_),
    psdProperties_(p.psdProperties_)
{
    initStatisticalSampling();
}

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "AerosolPopeParticleIO.C"

// ************************************************************************* //
