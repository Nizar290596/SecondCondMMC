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

#include "SootPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Variables * * * * * * * * * * * * //

template<class ParticleType>
const Foam::wordList Foam::SootPopeParticle<ParticleType>::sootContribsNames_
{
    "ySootInC2H2",
    "ySootSgC2H2",
    "ySootOxO2",
    "ySootOxOH",
    "nSootInC2H2",
    "nSootCg"
};

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::SootPopeParticle<ParticleType>::setCellValues
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
void Foam::SootPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    // First update the thermodynamic state
    ParticleType::calc(cloud, td, dt, cellI);

    // Update the soot variables: ySoot, nSoot, sootVf, sootIm, WL
    cloud.soot().calSoot
    (
        this->Y(),
        this->T(),
        this->ySoot(),
        this->nSoot(),
        this->sootVf(),
        this->sootIm(),
        cloud.composition().particleMixture(this->Y()).rho(this->pc(), this->T()),
        this->sootContribs()
    );
    
    //- Update the XO2XN2
    cloud.calcXO2XN2(this->Y(), this->XO2XN2());

    fVsoot() = 0.;

    Intermittency() = 1.;

    // Update soot
    //fVsoot() = cloud.calcfVsoot(this->Y(),this->T(),this->pc());

    //Intermittency() = cloud.calcIntermittency(fVsoot());


    // TODO: 1. add a function calculating the soot volume fraction
    //       2. add a function for soot intermittency
    //       3. then, update the radiation and T based on sootVf_
    // ---------------------------------Notes-------------------------------
    // it can be seen that, above functions are called from the cloud class,
    // thus they need to be defined in cloud class or submodel class,
    // here a new sootModel subclass has been built, it will be implemented
    // in the submodel. WL
    // ---------------------------------------------------------------------

    // Update enthalpy based on radiation using first fVsoot data
    this->hA() += cloud.radiation().Qrad
                    (
                        this->Y(),this->sootVf(),this->pc(),this->T()
                    )*dt;

    // update the temperature 
    ParticleType::calc(cloud, td, dt, cellI);
}


template<class ParticleType>
void Foam::SootPopeParticle<ParticleType>::mixProperties
(
    SootPopeParticle<ParticleType>& p, 
    SootPopeParticle<ParticleType>& q, 
    scalar mixExtent
)
{
    ParticleType::mixProperties(p,q,mixExtent);
    
}


template<class ParticleType>
void Foam::SootPopeParticle<ParticleType>::mixProperties
(
    SootPopeParticle<ParticleType>& p, 
    SootPopeParticle<ParticleType>& q, 
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{
    // First mix properties without soot
    ParticleType::mixProperties(p,q,mixExtent);
    
    // Mix soot properties
    
    // currently, only soot number density and mass fraction are considered,
    // which is based on the two-equation soot model
    
    // soot number density
    scalar nSootAv = (p.wt() * p.nSoot() + q.wt() * q.nSoot())/(p.wt() + q.wt());
    p.nSoot() = p.nSoot() + mixExtentSoot * (nSootAv - p.nSoot());
    q.nSoot() = q.nSoot() + mixExtentSoot * (nSootAv - q.nSoot());

    // soot mass fraction
    scalar ySootAv = (p.wt() * p.ySoot() + q.wt() * q.ySoot())/(p.wt() + q.wt());
    p.ySoot() = p.ySoot() + mixExtentSoot * (ySootAv - p.ySoot());
    q.ySoot() = q.ySoot() + mixExtentSoot * (ySootAv - q.ySoot());
}


template<class ParticleType>
void Foam::SootPopeParticle<ParticleType>::mixProperties
(
    SootPopeParticle<ParticleType>& p, 
    SootPopeParticle<ParticleType>& q, 
    scalar mixExtent,
    scalarList ScaledExtent
)
{
    ParticleType::mixProperties(p,q,mixExtent,ScaledExtent);
}


template<class ParticleType>
template<class CloudType>
void Foam::SootPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::SootPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::SootPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::SootPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();
    
    this->nameVariableLookUpTable().addNamedVariable("nSoot",nSoot_);
    this->nameVariableLookUpTable().addNamedVariable("sootVf",sootVf_);
    this->nameVariableLookUpTable().addNamedVariable("ySoot",ySoot_);

    // Add reference variables
    forAll(sootContribsNames_,i)
    {
        this->nameVariableLookUpTable().addNamedVariable
        (
            sootContribsNames_[i],
            sootContribs_[i]
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::SootPopeParticle<ParticleType>::SootPopeParticle
(
    const SootPopeParticle<ParticleType>& p
)
    :
    ParticleType(p),
    fVsoot_(p.fVsoot_),
    Intermittency_(p.Intermittency_),
    sootVf_(p.sootVf_),
    sootIm_(p.sootIm_),
    ySoot_(p.ySoot_),
    nSoot_(p.nSoot_),
    sootContribs_(p.sootContribs_),
    XO2XN2_(p.XO2XN2_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::SootPopeParticle<ParticleType>::SootPopeParticle
(
    const SootPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
    :
    ParticleType(p, mesh),
    fVsoot_(p.fVsoot_),
    Intermittency_(p.Intermittency_),
    sootVf_(p.sootVf_),
    sootIm_(p.sootIm_),
    ySoot_(p.ySoot_),
    nSoot_(p.nSoot_),
    sootContribs_(p.sootContribs_),
    XO2XN2_(p.XO2XN2_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SootPopeParticleIO.C"

// ************************************************************************* //
