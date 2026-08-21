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
inline Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
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
    secondCondFlag_(0),
    omegaOU_(0.0),
    phi_(0.0),
    phiModified_(0.0)
{}


// * * * * * * * * * * MixingPopeParticle Member Functions * * * * * * * * * * * * //

template<class ParticleType>
inline const Foam::scalarField& Foam::MixingPopeParticle<ParticleType>::XiR() const
{
    return XiR_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::MixingPopeParticle<ParticleType>::XiR
(
    const word& name
) const//!!!
{
    return XiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline const Foam::scalarField& Foam::MixingPopeParticle<ParticleType>::dXiR() const
{
    return dXiR_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::MixingPopeParticle<ParticleType>::dXiR
(
    const word& name
) const//!!!
{
    return dXiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalar Foam::MixingPopeParticle<ParticleType>::dx() const
{
    return dx_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::MixingPopeParticle<ParticleType>::XiR()
{
    return XiR_;
}


template<class ParticleType>
inline Foam::scalar& Foam::MixingPopeParticle<ParticleType>::XiR
(
    const word& name
)//!!!
{
    return XiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalarField& Foam::MixingPopeParticle<ParticleType>::dXiR()
{
    return dXiR_;
}


template<class ParticleType>
inline Foam::scalar& Foam::MixingPopeParticle<ParticleType>::dXiR
(
    const word& name
)//!!!
{
    return dXiR_[indexInXiR_[name]];
}


template<class ParticleType>
inline Foam::scalar& Foam::MixingPopeParticle<ParticleType>::dx()
{
    return dx_;
}

// Second conditioning accessors
 
template<class ParticleType>
inline Foam::label
Foam::MixingPopeParticle<ParticleType>::secondCondFlag() const
{
    return secondCondFlag_;
}
 
 
template<class ParticleType>
inline Foam::label&
Foam::MixingPopeParticle<ParticleType>::secondCondFlag()
{
    return secondCondFlag_;
}
 
 
template<class ParticleType>
inline Foam::scalar
Foam::MixingPopeParticle<ParticleType>::omegaOU() const
{
    return omegaOU_;
}
 
 
template<class ParticleType>
inline Foam::scalar&
Foam::MixingPopeParticle<ParticleType>::omegaOU()
{
    return omegaOU_;
}
 
template<class ParticleType>
inline Foam::scalar
Foam::MixingPopeParticle<ParticleType>::phi() const
{
    return phi_;
}
 
 
template<class ParticleType>
inline Foam::scalar&
Foam::MixingPopeParticle<ParticleType>::phi()
{
    return phi_;
} 

template<class ParticleType>
inline Foam::scalar
Foam::MixingPopeParticle<ParticleType>::phiModified() const
{
    return phiModified_;
}
 
 
template<class ParticleType>
inline Foam::scalar&
Foam::MixingPopeParticle<ParticleType>::phiModified()
{
    return phiModified_;
}

// ************************************************************************* //

