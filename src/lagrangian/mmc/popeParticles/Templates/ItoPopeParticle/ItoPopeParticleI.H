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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class ParticleType>
inline Foam::ItoPopeParticle<ParticleType>::ItoPopeParticle
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:

    ParticleType(owner, coordinates, cellI, tetFaceI, tetPtI),

    sCell_(),
    
    mass_(0),
    
    dpMsource_(0),
    
    wt_(0),
    
    A_(vector::zero),
    
    rw_(vector::zero),
    
    toPos_(this->position())
{}


template<class ParticleType>
inline Foam::ItoPopeParticle<ParticleType>::ItoPopeParticle
(
    const polyMesh& owner,
//    const vector& position,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const scalar& mass,
    const scalar& dpMsource,
    const scalar& wt,
    const vector &A
)
:
    ParticleType(owner, coordinates, cellI, tetFaceI, tetPtI), 
    
    sCell_(),
    
    mass_(mass),
    
    dpMsource_(dpMsource),
    
    wt_(wt),
    
    A_(A),
    
    rw_(vector::zero),
    
    toPos_(this->position())

{}


template<class ParticleType>
inline void Foam::ItoPopeParticle<ParticleType>::updateDispersedSources()
{
    //- Add interphase source terms (zero for gaseous combustion)
    this->m() += dpMsource();
    
    // Update of weights could be done outside of particle or deltaM needs to be 
    // included as a static member of the itoPopeParticle
        
    // this->wt() = this->m();
         
    dpMsource() = 0;
}    
// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //


// ************************************************************************* //
