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

template<class ParticleType>
inline Foam::PremixedMixingPopeParticle<ParticleType>::PremixedMixingPopeParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParticleType(mesh, coordinates, cellI, tetFaceI, tetPtI),
    XiR_(0),
    dXiR_(0),
    dx_(0.0),
    mixTime_(0.0),
    mixExt_(0.0),
    progVarInt_(0.0),
    XiRrlxFactor_(0.0),
    rlxStep_(0)
{
    initStatisticalSampling();
}


// * * * * * * * * * * PremixedMixingPopeParticle Member Functions * * * * * * * * * * * * //

template<class ParticleType>
inline const Foam::scalarField& Foam::PremixedMixingPopeParticle<ParticleType>::XiR() const
{
    return XiR_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::XiR
(
    const word& name
) const//!!!
{
    return XiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline const Foam::scalarField& Foam::PremixedMixingPopeParticle<ParticleType>::dXiR() const
{
    return dXiR_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::dXiR
(
    const word& name
) const//!!!
{
    return dXiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::dx() const
{
    return dx_;
}


template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::mixTime() const
{
    return mixTime_;
}


template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::mixExt() const
{
    return mixExt_;
}

template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::progVarInt() const
{
    return progVarInt_;
}


template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::XiRrlxFactor() const
{
    return XiRrlxFactor_;
}


template<class ParticleType>
inline Foam::scalar Foam::PremixedMixingPopeParticle<ParticleType>::rlxStep() const
{
    return rlxStep_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::PremixedMixingPopeParticle<ParticleType>::XiR()
{
    return XiR_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::XiR
(
    const word& name
)//!!!
{
    return XiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalarField& Foam::PremixedMixingPopeParticle<ParticleType>::dXiR()
{
    return dXiR_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::dXiR
(
    const word& name
)//!!!
{
    return dXiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::dx()
{
    return dx_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::mixTime()
{
    return mixTime_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::mixExt()
{
    return mixExt_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::progVarInt()
{
    return progVarInt_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::XiRrlxFactor()
{
    return XiRrlxFactor_;
}


template<class ParticleType>
inline Foam::scalar& Foam::PremixedMixingPopeParticle<ParticleType>::rlxStep()
{
    return rlxStep_;
}

// ************************************************************************* //
