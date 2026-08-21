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

#include "ThermoPopeParticle.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::ThermoPopeParticle<ParticleType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::setCellValues(cloud, td, dt, cellI);
    
    tetIndices tetIs = this->currentTetIndices();
    
    pc() = td.pInterp().interpolate(this->coordinates(),tetIs);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::ThermoPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    ParticleType::calc(cloud, td, dt, cellI);

    T() = cloud.composition().particleMixture(this->Y()).THa
                    (
                        this->hA(), this->pc(), this->T()
                    );

    hEqv() = cloud.composition().particleMixture(this->Y()).Hs
                    (
                        this->pc(),this->T()
                    );
}


template<class ParticleType>
void Foam::ThermoPopeParticle<ParticleType>::mixProperties
(
    ThermoPopeParticle<ParticleType>& p, 
    ThermoPopeParticle<ParticleType>& q, 
    scalar mixExtent
)
{
//    ParticleType::mixProperties(p,q,mixExtent);
    
    // Mix coupling variables

        forAll(p.XiC(),nn)
        { 
            //- Weighted pair mean
            scalar XiCAv = (p.wt() * p.XiC()[nn] + q.wt() * q.XiC()[nn])
                           /(p.wt() + q.wt());

            //- Mix them
            p.XiC()[nn] = p.XiC()[nn] + mixExtent * (XiCAv - p.XiC()[nn]);
            q.XiC()[nn] = q.XiC()[nn] + mixExtent * (XiCAv - q.XiC()[nn]);
        }

    // Mix enthalpy

        scalar hAv = (p.wt() * p.hA() + q.wt() * q.hA())/(p.wt() + q.wt());
            
        p.hA() = p.hA() + mixExtent * (hAv - p.hA());
        q.hA() = q.hA() + mixExtent * (hAv - q.hA());

    // Mix species
           
        //- Weighted pair mean
        scalarField YAv = (p.wt() * p.Y() + q.wt() * q.Y())/(p.wt() + q.wt());

        //- Mix them
            p.Y() = p.Y() + mixExtent * (YAv - p.Y());
            q.Y() = q.Y() + mixExtent * (YAv - q.Y());
}

template<class ParticleType>
void Foam::ThermoPopeParticle<ParticleType>::mixProperties
(
    ThermoPopeParticle<ParticleType>& p, 
    ThermoPopeParticle<ParticleType>& q, 
    scalar mixExtent,
    scalarList ScaledExtent
)
{
//    ParticleType::mixProperties(p,q,mixExtent);
    
    // Mix coupling variables

        forAll(p.XiC(),nn)
        { 
            //- Weighted pair mean
            scalar XiCAv = (p.wt() * p.XiC()[nn] + q.wt() * q.XiC()[nn])
                           /(p.wt() + q.wt());

            //- Mix them
            p.XiC()[nn] = p.XiC()[nn] + mixExtent * (XiCAv - p.XiC()[nn]);
            q.XiC()[nn] = q.XiC()[nn] + mixExtent * (XiCAv - q.XiC()[nn]);
        }

    // Mix enthalpy

        scalar hAv = (p.wt() * p.hA() + q.wt() * q.hA())/(p.wt() + q.wt());
            
        p.hA() = p.hA() + mixExtent * (hAv - p.hA());
        q.hA() = q.hA() + mixExtent * (hAv - q.hA());

    // Mix species
           
        //- Weighted pair mean
        scalarField YAv = (p.wt() * p.Y() + q.wt() * q.Y())/(p.wt() + q.wt());

        scalar sumYp(0);
        scalar sumYq(0);


            //Pout<< "p.Y() before mixing before scaled=" << p.Y() << endl;
            //Pout<< "sumYp before mixing before scaled=" << sumYp << endl;

        const scalarList b4Mix_pY = p.Y();
        const scalarList b4Mix_qY = q.Y();

        //- Mix them
        forAll(p.Y(),Spi)
        {
            p.Y()[Spi] = p.Y()[Spi] + ScaledExtent[Spi] * (YAv[Spi] - p.Y()[Spi]);
            q.Y()[Spi] = q.Y()[Spi] + ScaledExtent[Spi] * (YAv[Spi] - q.Y()[Spi]);

            sumYp += p.Y()[Spi];
            sumYq += q.Y()[Spi];
        }

        const scalarList afterMixb4scaled_pY = p.Y();
        const scalarList afterMixb4scaled_qY = q.Y();


        //normalise Y
        scalar norFactor_p = 1.0/ sumYp;
        scalar norFactor_q = 1.0/ sumYq;

        //sumYp = 0;
        //sumYq = 0;

        forAll(p.Y(),Spi)
        {
            p.Y()[Spi] = norFactor_p * p.Y()[Spi];
            q.Y()[Spi] = norFactor_q * q.Y()[Spi];

           // sumYp += p.Y()[Spi];
           // sumYq += q.Y()[Spi];

        }

    if (max(ScaledExtent) > 1.0)
    {
            Pout << "ScaledExtent ="<< ScaledExtent << endl;
        FatalErrorIn
                    (
                        "  "
                    )   
                        << "max(ScaledExtent) > 1.0" << nl
                        << exit(FatalError);
    }



    //Fatal error if the scaled mixing cause 1% mass overshoot.
    if(sumYp <0.99 ||sumYp>1.01 || sumYq<0.99 || sumYq>1.01)
    {
        Pout << "p.Y() before mixing =" << b4Mix_pY << endl;
        Pout << "p.Y() after mixing before scaled=" << afterMixb4scaled_pY << endl;
        Pout << "p.Y() after mixing after scaled=" << p.Y() << endl;
        Pout << "sumYp after mixing before scaled=" << sumYp << endl;
        Pout << "norFactor_p ="<< norFactor_p << endl;

        Pout << "q.Y() before mixing =" << b4Mix_qY << endl;
        Pout << "q.Y() after mixing before scaled=" << afterMixb4scaled_qY << endl;
        Pout << "q.Y() after mixing after scaled=" << q.Y() << endl;
        Pout << "sumYq after mixing before scaled=" << sumYq << endl;
        Pout << "norFactor_q ="<< norFactor_q << endl;

        Pout << "ScaledExtent ="<< ScaledExtent << endl;

    FatalErrorIn
                (
                    "  "
                )   
                    << "crash at mixing" << nl
                    << exit(FatalError);

    }
}


template<class ParticleType>
template<class CloudType>
void Foam::ThermoPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
    ParticleType::setStaticProperties(c);
    
    indexInXiC_ = c.coupling().XiC().cVarInXiC();
    XiCNames_ = c.coupling().XiCNames();
    componentNames_ = c.composition().componentNames();
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::ThermoPopeParticle<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalData(vars);
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::ThermoPopeParticle<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    return ParticleType::getStatisticalDataNames(vars);
}


template<class ParticleType>
void Foam::ThermoPopeParticle<ParticleType>::initStatisticalSampling()
{
    ParticleType::initStatisticalSampling();

    // Add ethalpy
    this->nameVariableLookUpTable().addNamedVariable("hEqv",hEqv_);
    this->nameVariableLookUpTable().addNamedVariable("T",T_);

    // Note: XiCNames is initialized as a static variablebefore 
    // XiC_ is set. Therefore we need to check if they are set to avoid
    // accessing elements out of bounds

    if (XiCNames_.size() == XiC_.size())
    {
        // Add reference variables
        for (const word& name : XiCNames_)
        {
            this->nameVariableLookUpTable().addNamedVariable
            (
                name,
                XiC(name)
            );
        }
    }

    // Add species
    if (componentNames_.size() == Y_.size())
    {
        forAll(componentNames_,i)
        {
            this->nameVariableLookUpTable().addNamedVariable(componentNames_[i],Y_[i]);
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ThermoPopeParticle<ParticleType>::ThermoPopeParticle
(
    const ThermoPopeParticle<ParticleType>& p
)
    :
    ParticleType(p),
    T_(p.T_),
    Y_(p.Y_),
    XiC_(p.XiC_),
    hA_(p.hA_),
    hEqv_(p.hEqv_),
    dpYsource_(p.dpYsource_),
    dpXiCsource_(p.dpXiCsource_),
    dp_hAsource_(p.dp_hAsource_),
    pc_(p.pc_),
    NumActSp_(p.NumActSp_)
{
    initStatisticalSampling();
}


template<class ParticleType>
Foam::ThermoPopeParticle<ParticleType>::ThermoPopeParticle
(
    const ThermoPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
    :
    ParticleType(p, mesh),
    T_(p.T_),
    Y_(p.Y_),
    XiC_(p.XiC_),
    hA_(p.hA_),
    hEqv_(p.hEqv_),
    dpYsource_(p.dpYsource_),
    dpXiCsource_(p.dpXiCsource_),
    dp_hAsource_(p.dp_hAsource_),
    pc_(p.pc_),
    NumActSp_(p.NumActSp_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ThermoPopeParticleIO.C"

// ************************************************************************* //
