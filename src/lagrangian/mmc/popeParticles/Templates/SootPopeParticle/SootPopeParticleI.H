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
inline Foam::SootPopeParticle<ParticleType>::SootPopeParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParticleType(mesh, coordinates, cellI, tetFaceI, tetPtI),
    fVsoot_(0),
    Intermittency_(0),
    sootVf_(0),
    sootIm_(0),
    ySoot_(0),
    nSoot_(0),
    sootContribs_(6,0.0),
    XO2XN2_(0)
{}



// * * * * * * * * * * SootPopeParticle Member Functions * * * * * * * * * * * * //

template<class ParticleType>
inline const Foam::scalarField& Foam::SootPopeParticle<ParticleType>::fVsoot() const
{
    return fVsoot_;
}


template<class ParticleType>
inline const Foam::scalarField& Foam::SootPopeParticle<ParticleType>::Intermittency() const
{
    return Intermittency_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::SootPopeParticle<ParticleType>::sootVf() const
{
    return sootVf_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::SootPopeParticle<ParticleType>::sootIm() const
{
    return sootIm_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::SootPopeParticle<ParticleType>::ySoot() const
{
    return ySoot_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::SootPopeParticle<ParticleType>::nSoot() const
{
    return nSoot_;
}


template<class ParticleType>
inline const Foam::scalarList& Foam::SootPopeParticle<ParticleType>::sootContribs() const
{
    return sootContribs_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::SootPopeParticle<ParticleType>::XO2XN2() const
{
    return XO2XN2_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::SootPopeParticle<ParticleType>::fVsoot()
{
    return fVsoot_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::SootPopeParticle<ParticleType>::Intermittency()
{
    return Intermittency_;
}


template<class ParticleType>
inline  Foam::scalar& Foam::SootPopeParticle<ParticleType>::sootVf() 
{
    return sootVf_;
}


template<class ParticleType>
inline Foam::scalar& Foam::SootPopeParticle<ParticleType>::sootIm() 
{
    return sootIm_;
}


template<class ParticleType>
inline Foam::scalar& Foam::SootPopeParticle<ParticleType>::ySoot() 
{
    return ySoot_;
}


template<class ParticleType>
inline Foam::scalar& Foam::SootPopeParticle<ParticleType>::nSoot() 
{
    return nSoot_;
}


template<class ParticleType>
inline Foam::scalarList& Foam::SootPopeParticle<ParticleType>::sootContribs() 
{
    return sootContribs_;
}


template<class ParticleType>
inline Foam::scalar& Foam::SootPopeParticle<ParticleType>::XO2XN2() 
{
    return XO2XN2_;
}

template<class ParticleType>
inline void Foam::SootPopeParticle<ParticleType>::updateDispersedSources()
{
    ParticleType::updateDispersedSources();
}


// ************************************************************************* //
