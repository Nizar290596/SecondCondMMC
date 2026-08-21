/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    T0_(this->dict_, 0.0),
    TMin_(this->dict_, 0.0),
    TMax_(this->dict_, VGREAT),
    Cp0_(this->dict_, 0.0),
    epsilon0_(this->dict_, 0.0),
    f0_(this->dict_, 0.0),
    Tvap0_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    TMax_(cp.TMax_),
    Cp0_(cp.Cp0_),
    epsilon0_(cp.epsilon0_),
    f0_(cp.f0_),
    Tvap0_(cp.Tvap0_)
{}


template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    T0_(this->dict_, "T0"),
    TMin_(this->dict_, "TMin", 200.0),
    TMax_(this->dict_, "TMax", 5000.0),
    Cp0_(this->dict_, "Cp0"),
    epsilon0_(this->dict_, "epsilon0"),
    f0_(this->dict_, "f0"),
    Tvap0_(this->dict_,"Tvap0")
{}


template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::DropletSprayDNSThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    T_(0.0),
    Cp_(0.0),
    mass0_(this->mass()),    
    jFuel_(-1),
    Tvap_(0.0)
{
    initStatisticalSampling();
}


template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::DropletSprayDNSThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    T_(0.0),
    Cp_(0.0),
    mass0_(this->mass()),    
    jFuel_(-1),
    Tvap_(0.0)
{
    initStatisticalSampling();
}


template<class ParcelType>
inline Foam::DropletSprayDNSThermoParcel<ParcelType>::DropletSprayDNSThermoParcel
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
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
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
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    T_(constProps.T0()),
    Cp_(constProps.Cp0()),
    mass0_(this->mass()),   
    jFuel_(-1),
    Tvap_(constProps.Tvap0())
{
    initStatisticalSampling();
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //


template<class ParcelType>
inline void Foam::DropletSprayDNSThermoParcel<ParcelType>::reset()
{
    mDot_ = 0;
    ParcelType::reset();
}

template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::TMax() const
{
    return TMax_.value();
}


template<class ParcelType>
inline void
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::setTMax(const scalar TMax)
{
    TMax_.setValue(TMax);
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::Cp0() const
{
    return Cp0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::Tvap0() const
{
    return Tvap0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::DropletSprayDNSThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_.value();
}

// * * * * * * * * * * DropletSprayDNSThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}

template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::mass0() const
{
    return mass0_;
}


template<class ParcelType>
inline Foam::label Foam::DropletSprayDNSThermoParcel<ParcelType>::jFuel() const
{
    return jFuel_;
}


template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::hs() const
{
    return Cp_*(T_ - 298.15);
}

template<class ParcelType>
inline Foam::scalar Foam::DropletSprayDNSThermoParcel<ParcelType>::Tvap() const
{
    return Tvap_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::T()
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::Cp()
{
    return Cp_;
}

template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::Tvap()
{
    return Tvap_;
}

template<class ParcelType>
inline Foam::scalar& Foam::DropletSprayDNSThermoParcel<ParcelType>::mass0()
{
    return mass0_;
}

// ************************************************************************* //
