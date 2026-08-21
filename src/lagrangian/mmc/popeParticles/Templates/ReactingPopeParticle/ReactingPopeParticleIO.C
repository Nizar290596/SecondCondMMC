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

#include "ReactingPopeParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReactingPopeParticle<ParticleType>::ReactingPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
: 
    ParticleType(mesh, is, readFields, newFormat)
{
    if (readFields)
        is >> cpuTime_;
}


template<class ParticleType>
template<class CloudType>
void Foam::ReactingPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
      return;
    }

    ParticleType::readFields(c);
}


template<class ParticleType>
template<class CloudType>
void Foam::ReactingPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.cpuTime();
    
    // Check state of Ostream
    os.check
    (
        "operator<<(Ostream& os, const popeParticle&)"
    );

    return os;
}


// ************************************************************************* //
