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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//using namespace Foam;

template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    T0_(this->dict_, 0.0),
    TMin_(this->dict_, 0.0),
    TMax_(this->dict_, VGREAT),
    Sc_(this->dict_, 0.0),
    Le_(this->dict_, 0.0),
    LH_(this->dict_, 0.0),
    FH_(this->dict_, 0.0),
    molWtSstate0_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    TMax_(cp.TMax_),
    Sc_(cp.Sc_),
    Le_(cp.Le_),
    LH_(cp.LH_),
    FH_(cp.FH_),
    molWtSstate0_(cp.molWtSstate0_)
{}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    ParcelType::constantProperties(parentDict),
    T0_(this->dict_, "T0"),
    TMin_(this->dict_, "TMin"),
    TMax_(this->dict_, "TMax"),
    Sc_(this->dict_, "Sc"),
    Le_(this->dict_, "Le"),
    LH_(this->dict_, "latentHeat"),
    FH_(this->dict_, "formationHeat", 0.0),
    molWtSstate0_(this->dict_, "molWtSstate0")
{}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar minParticleMass,
    const scalar T0,
    const scalar TMin,
    const scalar TMax,
    const scalar Sc,
    const scalar Le,
    const scalar latentHeat,
    const scalar formationHeat,
    const scalar molWtSstate0
)
:
    ParcelType::constantProperties
    (
        parcelTypeId,
        rhoMin,
        minParticleMass
    ),
    T0_(T0),
    TMin_(TMin),
    TMax_(TMax),
    Sc_(Sc),
    Le_(Le),
    LH_(latentHeat),
    FH_(formationHeat),
    molWtSstate0_(molWtSstate0)
{}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    TD_(0.0),
    Cp_(0.0),
    TGas_(0.0),
    TFilm_(0.0),
    TSurf_(0.0),
    pc_(0.0),
    fGas_(0.0),
    fSurf_(0.0),
    fFilm_(0.0),
    YFuel_(3),
    YFuelVap_(3),
    YGas_(3),
    YFilm_(3),
    YSurf_(3),
    haGas_(0.0),
    hD_(0.0),
    haSurf_(0.0),
    haFilm_(0.0),
    haFuelVap_(0.0),
    muFilm_(0.0),
    rhoFilm_(0.0),
    qR_(0.0),
    qD_(0.0),
    mFlux_(0.0),
    molWtSstate_(0.0),   
    ZFilm_(0.0),
    BM_(0.0),
    superCelli_(0.0)
{
    initStatisticalSampling();
}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    TD_(0.0),
    Cp_(0.0),
    TGas_(0.0),
    TFilm_(0.0),
    TSurf_(0.0),
    pc_(0.0),
    fGas_(0.0),
    fSurf_(0.0),
    fFilm_(0.0),
    YFuel_(3),
    YFuelVap_(3),
    YGas_(3),
    YFilm_(3),
    YSurf_(3),
    haGas_(0.0),
    hD_(0.0),
    haSurf_(0.0),
    haFilm_(0.0),
    haFuelVap_(0.0),
    muFilm_(0.0),
    rhoFilm_(0.0),
    qR_(0.0),
    qD_(0.0),
    mFlux_(0.0),
    molWtSstate_(0.0),   
    ZFilm_(0.0),
    BM_(0.0),
    superCelli_(0.0)
{
    initStatisticalSampling();
}


template<class ParcelType>
inline Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        constProps
    ),
    TD_(constProps.T0()),
    Cp_(0.0),
    TGas_(0.0),
    TFilm_(0.0),
    TSurf_(0.0),
    pc_(0.0),
    fGas_(0.0),
    fSurf_(0.0),
    fFilm_(0.0),
    YFuel_(3),
    YFuelVap_(3),
    YGas_(3),
    YFilm_(3),
    YSurf_(3),
    haGas_(0.0),
    hD_(0.0),
    haSurf_(0.0),
    haFilm_(0.0),
    haFuelVap_(0.0),
    muFilm_(0.0),
    rhoFilm_(0.0),
    qR_(0.0),
    qD_(0.0),
    mFlux_(0.0),
    molWtSstate_(constProps.molWtSstate0()),   
    ZFilm_(0.0),
    BM_(0.0),
    superCelli_(0.0)
{
    initStatisticalSampling();
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::TMax() const
{
    return TMax_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::Sc() const
{
    return Sc_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::Le() const
{
    return Le_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::LH() const
{
    return LH_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::FH() const
{
    return FH_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::sprayThermoParcel<ParcelType>::constantProperties::molWtSstate0() const
{
    return molWtSstate0_.value();
}

// * * * * * * * * * * sprayThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::T() const
{
    return TD_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::hs() const
{
    // Equation taken from thermoParcel of OpenFOAM v2206
    return Cp_*(T() - 298.15);
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::TGas() const
{
    return TGas_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::TFilm() const
{
    return TFilm_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::TSurf() const
{
    return TSurf_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::hD() const
{
    return hD_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::pc() const
{
    return pc_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::fGas() const
{
    return fGas_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::fSurf() const
{
    return fSurf_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::fFilm() const
{
    return fFilm_;
}


template<class ParcelType>
inline Foam::label Foam::sprayThermoParcel<ParcelType>::nThermoSpecies() const
{
    return nThermoSpecies_;
}


template<class ParcelType>
inline Foam::label Foam::sprayThermoParcel<ParcelType>::nFuelSpecies() const
{
    return nFuelSpecies_;
}


template<class ParcelType>
inline Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YFuel() const
{
    return YFuel_;
}


template<class ParcelType>
inline Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YFuelVap() const
{
    return YFuelVap_;
}


template<class ParcelType>
inline Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YFilm() const
{
    return YFilm_;
}


template<class ParcelType>
inline Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YSurf() const
{
    return YSurf_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::haGas() const
{
    return haGas_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::haSurf() const
{
    return haSurf_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::haFilm() const
{
    return haFilm_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::haFuelVap() const
{
    return haFuelVap_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::rhoFilm() const 
{
    return rhoFilm_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::muFilm() const
{
    return muFilm_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::qR() const
{
    return qR_;
}

template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::qD() const
{
    return qD_;
}

template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::mFlux() const
{
    return mFlux_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::molWtSstate() const
{
    return molWtSstate_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::ZFilm() const
{
    return ZFilm_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::BM() const
{
    return BM_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::PrFilm
(   
    const scalar muFilm,
    const scalar rhoFilm,
    const scalar alpha
) const    
{
    scalar PrandtlNumber = (muFilm/rhoFilm)/alpha;
    return PrandtlNumber;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::NuFilm
(   
    const scalar ReFilm,
    const scalar PrFilm
) const 
{
    scalar NusseltNumber = 2.0 + 0.6*sqrt(ReFilm)*pow(PrFilm,1.0/3.0);
    return NusseltNumber;
}

    
template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::ShFilm
(   
    const scalar ReFilm,
    const scalar Sc
) const
{
    scalar SherwoodNumber = 2.0 + 0.6*sqrt(ReFilm)*pow(Sc, 1.0/3.0);
    return SherwoodNumber;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::BM 
(
    const scalar fGas,           // Gas mixture fraction
    const scalar fSurf           // Surface surface mixture fraction
) const
{
    scalar bmNumber = ((fGas-fSurf)/(fSurf -1.0+SMALL));
    return bmNumber;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::BH
(
    const scalar haGas,           // Gas enthalpy
    const scalar haSurf,          // Surface enthalpy
    const scalar hD,              // Dispersion phase enthalpy 
    const scalar qD,              // Internal heating  
    const scalar qR               // heat due to radiation   
)const
{
    scalar bhNumber = ((haGas -haSurf)/(haSurf - hD + qD -qR)); 
    return bhNumber;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::qR
(
    const scalar TD,          // dispersion phase temperature
    const scalar TWall,       // gas temperature
    const scalar mFlux        // mass flux
) const
{
    scalar qRadiation = -5.6703E-08*(pow(TD, 4.0) - pow(TWall, 4.0))/fabs(mFlux);
    return qRadiation;
}


//- edit functions
template<class ParcelType>
inline Foam::scalar& Foam::sprayThermoParcel<ParcelType>::T()
{
    return TD_;
}

template<class ParcelType>
inline Foam::scalar& Foam::sprayThermoParcel<ParcelType>::Cp()
{
     return Cp_;
}

template<class ParcelType>
inline Foam::scalar& Foam::sprayThermoParcel<ParcelType>::TGas()
{
    return TGas_;
}


template<class ParcelType>
inline Foam::label& Foam::sprayThermoParcel<ParcelType>::nThermoSpecies()
{
    return nThermoSpecies_;
}


template<class ParcelType>
inline Foam::label& Foam::sprayThermoParcel<ParcelType>::nFuelSpecies()
{
    return nFuelSpecies_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::sprayThermoParcel<ParcelType>::YFuel()
{
    return YFuel_;
}


template<class ParcelType>
inline Foam::scalar& Foam::sprayThermoParcel<ParcelType>::molWtSstate() 
{
    return molWtSstate_;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::temperatureUpdateSwitch
(
    scalar F0, 
    scalar F1
)
{
    scalar ratio = 0.5; // ratio of Crank Nicolsen method 
    scalar dTdt = ratio * F0 + (1.0 - ratio) * F1;

    return dTdt;
}


template<class ParcelType>
inline Foam::scalar Foam::sprayThermoParcel<ParcelType>::calcMMCSource
(
    const scalar dt,
    const scalar c, 
    const scalar phi_n, 
    const scalar tau
)
{
    scalar dPhi = (1 - exp(-dt / tau)) * (c - phi_n);
    return dPhi;
}

// ************************************************************************* //
