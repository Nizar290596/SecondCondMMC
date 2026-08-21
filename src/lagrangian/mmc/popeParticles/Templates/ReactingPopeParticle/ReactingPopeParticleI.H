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
inline Foam::ReactingPopeParticle<ParticleType>::ReactingPopeParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParticleType(mesh, coordinates, cellI, tetFaceI, tetPtI)
{}

// * * * * * * * * * * ReactingPopeParticle Member Functions * * * * * * * * //

template<class ParticleType>
inline Foam::scalar& Foam::ReactingPopeParticle<ParticleType>::cpuTime()
{
    return cpuTime_;
}


template<class ParticleType>
inline const Foam::scalar& 
Foam::ReactingPopeParticle<ParticleType>::cpuTime() const
{
    return cpuTime_;
}
// ************************************************************************* //
