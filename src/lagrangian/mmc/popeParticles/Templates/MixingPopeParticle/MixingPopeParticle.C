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

#include "MixingPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::MixingPopeParticle<ParticleType>::setCellValues
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
              
            //if (setOfXi.Vars(indexInXi[*iter]).refType()=="interpolated")
            //    XiR()[indexInXiR[*iter]] = td.XiRInterp()[indexInXiR[*iter]].interpolate(this->coordinates(),tetIs);                
            //else if (setOfXi.Vars(indexInXi[*iter]).refType()=="evolved")
            //{
                // Do nothing since it is evolved by MMC model
		scalar XiROld = XiR(*iter);
		scalar XiRNew = setOfXi.Vars(indexInXi[*iter]).evolveMethod().compute(*this,dt,XiROld);
		XiR(*iter) = XiRNew;
		//Info << "XiROld " << XiROld << endl;
		//Info << "XiRNew " << XiRNew << endl;
           // }
           // else
           //     XiR()[indexInXiR[*iter]] = 0.0;//this->position();
        }
    }

}


template<class ParticleType>
template<class TrackCloudType>
void Foam::MixingPopeParticle<ParticleType>::calc
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
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    scalar mixExtent
)
{       
    //ParticleType::mixProperties(p, q, mixExtent);

//    const label fp = p.secondCondFlag();
//    const label fq = q.secondCondFlag();

    // S(φ): mix progress variable between the pair (same as other reactive
    // scalars in MMCcurl — weighted pair mean with mixExtent relaxation)
    //scalar phiAv = (p.wt() * p.phi() + q.wt() * q.phi()) / (p.wt() + q.wt());
    //p.phi() = p.phi() + mixExtent * (phiAv - p.phi());
    //q.phi() = q.phi() + mixExtent * (phiAv - q.phi());

    // Recompute φ° for flagged particles: phi has changed, ω_OU has not
    //    p.phiModified() = p.phi() * Foam::exp(secondCondBeta_s_ * p.omegaOU());
    //    q.phiModified() = q.phi() * Foam::exp(secondCondBeta_s_ * q.omegaOU());
//    if (fp == 0 && fq == 0)
//    {
        // Non-flagged pair: mix Y, T, hA via the parent chain.
        // φ is left untouched (evolves only via W(φ) for these particles).
//        ParticleType::mixProperties(p, q, mixExtent);
//    }
//    else if (fp == 1 && fq == 1)
//    {
        // Flagged pair: mix only φ here. Y, T, hA are mixed in step 5
        // by secondCondMMCcurl on the 4-D (φ°, sP) reference space, so
        // not mixing them here avoids double mixing of those scalars.
//        const scalar wtSum = p.wt() + q.wt();
//        if (wtSum > VSMALL)
//        {
//            const scalar phiAv = (p.wt()*p.phi() + q.wt()*q.phi())/wtSum;
//            p.phi() += mixExtent * (phiAv - p.phi());
//            q.phi() += mixExtent * (phiAv - q.phi());
            // phiModified is refreshed at step 4 (updateOUProcess).
//        }
//    }

     // Baseline first-conditioning composition mixing (mixes Y, hA, XiC via the
    // parent chain) for all pairs. When second conditioning is enabled,
    // MMCcurl::mixpair mixes only the progress variable phi instead and does
    // not call this function.
    ParticleType::mixProperties(p, q, mixExtent);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    const scalar& mixExtent,
    const scalar& mixExtentSoot
)
{       
    //ParticleType::mixProperties(p, q, mixExtent,mixExtentSoot);

    // S(φ): mix progress variable between the pair (same as other reactive
    // scalars in MMCcurl — weighted pair mean with mixExtent relaxation)
    //scalar phiAv = (p.wt() * p.phi() + q.wt() * q.phi()) / (p.wt() + q.wt());
    //p.phi() = p.phi() + mixExtent * (phiAv - p.phi());
    //q.phi() = q.phi() + mixExtent * (phiAv - q.phi());

    // Recompute φ° for flagged particles: phi has changed, ω_OU has not
    //    p.phiModified() = p.phi() * Foam::exp(secondCondBeta_s_ * p.omegaOU());
    //    q.phiModified() = q.phi() * Foam::exp(secondCondBeta_s_ * q.omegaOU());

//    const label fp = p.secondCondFlag();
//    const label fq = q.secondCondFlag();

//    if (fp == 0 && fq == 0)
//    {
//        ParticleType::mixProperties(p, q, mixExtent, mixExtentSoot);
//    }
//    else if (fp == 1 && fq == 1)
//    {
//        const scalar wtSum = p.wt() + q.wt();
//        if (wtSum > VSMALL)
//        {
//            const scalar phiAv = (p.wt()*p.phi() + q.wt()*q.phi())/wtSum;
//            p.phi() += mixExtent * (phiAv - p.phi());
//            q.phi() += mixExtent * (phiAv - q.phi());
//        }
//    }
    // Baseline first-conditioning composition mixing (soot-aware variant).
    // When second conditioning is enabled, MMCcurl::mixpair mixes only the
    // progress variable phi instead and does not call this function.
    ParticleType::mixProperties(p, q, mixExtent, mixExtentSoot);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::mixProperties
(
    MixingPopeParticle<ParticleType>& p,
    MixingPopeParticle<ParticleType>& q,
    scalar mixExtent,
    scalarList ScaledExtent
)
{
    //ParticleType::mixProperties(p, q, mixExtent,ScaledExtent);

    // S(φ): mix progress variable between the pair (same as other reactive
    // scalars in MMCcurl — weighted pair mean with mixExtent relaxation)
//    scalar phiAv = (p.wt() * p.phi() + q.wt() * q.phi()) / (p.wt() + q.wt());
//    p.phi() = p.phi() + mixExtent * (phiAv - p.phi());
//    q.phi() = q.phi() + mixExtent * (phiAv - q.phi());

    // Recompute φ° for flagged particles: phi has changed, ω_OU has not
//        p.phiModified() = p.phi() * Foam::exp(secondCondBeta_s_ * p.omegaOU());
//        q.phiModified() = q.phi() * Foam::exp(secondCondBeta_s_ * q.omegaOU());

//    const label fp = p.secondCondFlag();
//    const label fq = q.secondCondFlag();

//    if (fp == 0 && fq == 0)
//    {
//        ParticleType::mixProperties(p, q, mixExtent, ScaledExtent);
//    }
//    else if (fp == 1 && fq == 1)
//    {
//        const scalar wtSum = p.wt() + q.wt();
//        if (wtSum > VSMALL)
//        {
//            const scalar phiAv = (p.wt()*p.phi() + q.wt()*q.phi())/wtSum;
//            p.phi() += mixExtent * (phiAv - p.phi());
//            q.phi() += mixExtent * (phiAv - q.phi());
//        }
//    }
     // Baseline first-conditioning composition mixing (scaled-extent variant).
    // When second conditioning is enabled, MMCcurl::mixpair mixes only the
    // progress variable phi instead and does not call this function.
    ParticleType::mixProperties(p, q, mixExtent, ScaledExtent);
}


template<class ParticleType>
template<class CloudType>
void Foam::MixingPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);

    // Cache β so mixProperties (which has no cloud ref) can recompute φ°
    //secondCondBeta_s_ = c.secondCondBeta();

    if (c.mixing().numXiR())
    {
        indexInXiR_ = c.mixing().XiR().rVarInXiR();
        XiRNames_ = c.mixing().XiRNames();
    }
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::MixingPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::MixingPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::MixingPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();

    // Add ethalpy
    //this->nameVariableLookUpTable().addNamedVariable("dx",dx_);
    // Add physical-space mixing distance
    this->nameVariableLookUpTable().addNamedVariable("dx", dx_);

    // Second conditioning scalar
    this->nameVariableLookUpTable().addNamedVariable("omegaOU", omegaOU_);
    this->nameVariableLookUpTable().addNamedVariable("phi",	phi_);
    this->nameVariableLookUpTable().addNamedVariable("phiModified", phiModified_);

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
Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
(
    const MixingPopeParticle<ParticleType>& p
)
:
    ParticleType(p),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    secondCondFlag_(p.secondCondFlag_),
    omegaOU_(p.omegaOU_),
    phi_(p.phi_),
    phiModified_(p.phiModified_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
(
    const MixingPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
:
    ParticleType(p, mesh),
    XiR_(p.XiR_),
    dXiR_(p.dXiR_),
    dx_(p.dx_),
    secondCondFlag_(p.secondCondFlag_),
    omegaOU_(p.omegaOU_),
    phi_(p.phi_),
    phiModified_(p.phiModified_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MixingPopeParticleIO.C"


// ************************************************************************* //

