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
inline Foam::AerosolPopeParticle<ParticleType>::AerosolPopeParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParticleType(mesh, coordinates, cellI, tetFaceI, tetPtI),
    physicalAerosolProperties_(0),
    psdProperties_(0)
{}

// * * * * * * * * * * AerosolPopeParticle Member Functions * * * * * * * * * * * * //


template<class ParticleType>
inline const Foam::scalarField& 
Foam::AerosolPopeParticle<ParticleType>::physicalAerosolProperties() const
{
    return physicalAerosolProperties_;
}

template<class ParticleType>
inline const Foam::scalarField& 
Foam::AerosolPopeParticle<ParticleType>::psdProperties() const
{
    return psdProperties_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::AerosolPopeParticle<ParticleType>::physicalAerosolProperties()
{
    return physicalAerosolProperties_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::AerosolPopeParticle<ParticleType>::psdProperties()
{
    return psdProperties_;
}

// ************************************************************************* //
