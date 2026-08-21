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
inline Foam::ThermoPopeParticle<ParticleType>::ThermoPopeParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParticleType(mesh, coordinates, cellI, tetFaceI, tetPtI),
    T_(0.0),
    Y_(0),
    XiC_(0),
    hA_(0),
    hEqv_(0),
    dpYsource_(0),
    dpXiCsource_(0),
    dp_hAsource_(0),
    pc_(0.0),
    NumActSp_(0)
{}



// * * * * * * * * * * ThermoPopeParticle Member Functions * * * * * * * * * * * * //

template<class ParticleType>
inline Foam::scalar Foam::ThermoPopeParticle<ParticleType>::T() const
{
    return T_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::hA() const
{
    return hA_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::hEqv() const
{
    return hEqv_;
}


template<class ParticleType>
inline const Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::Y() const
{
    return Y_;
}


template<class ParticleType>
inline const Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::
dpYsource() const
{
    return dpYsource_;
}


template<class ParticleType>
inline const Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::
dpXiCsource() const
{
    return dpXiCsource_;
}


template<class ParticleType>
inline const Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::dpXiCsource
(
    const word& name
) const //!!!
{
    return dpXiCsource_[indexInXiC_[name]];
}


template<class ParticleType>
inline const Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::
dp_hAsource() const
{
    return dp_hAsource_;
}


template<class ParticleType>
inline const Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::
XiC() const //!!!
{
    return XiC_;
}

template<class ParticleType>
inline const Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::XiC
(
    const word& name
) const //!!!
{
    return XiC_[indexInXiC_[name]];
}

template<class ParticleType>
inline Foam::scalar Foam::ThermoPopeParticle<ParticleType>::pc() const
{
    return pc_;
}


template<class ParticleType>
inline const Foam::label& Foam::ThermoPopeParticle<ParticleType>::NumActSp() const
{
    return NumActSp_;
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::T()
{
    return T_;
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::hA()
{
    return hA_;
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::hEqv()
{
    return hEqv_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::Y()
{
    return Y_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::dpYsource()
{
    return dpYsource_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::dpXiCsource()
{
    return dpXiCsource_;
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::dpXiCsource
(
    const word& name
)//!!!
{
    return dpXiCsource_[indexInXiC_[name]];
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::dp_hAsource()
{
    return dp_hAsource_;
}


template<class ParticleType>
inline Foam::scalarField& Foam::ThermoPopeParticle<ParticleType>::XiC() //!!!
{
    return XiC_;
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::XiC
(
    const word& name
)//!!!
{
    return XiC_[indexInXiC_[name]];
}


template<class ParticleType>
inline Foam::scalar& Foam::ThermoPopeParticle<ParticleType>::pc() 
{
    return pc_;
}


template<class ParticleType>
inline Foam::label& Foam::ThermoPopeParticle<ParticleType>::NumActSp()
{
    return NumActSp_;
}


template<class ParticleType>
inline void Foam::ThermoPopeParticle<ParticleType>::updateDispersedSources()
{
    ParticleType::updateDispersedSources();

    //- Add interphase source terms (zero for gaseous combustion)
    this->Y() += dpYsource();
    
    this->XiC() += dpXiCsource();
    
    hA() += dp_hAsource();
    
    // Update of Temperature and density done outside in multiphase
    
    //- Reset interphase source terms to zero
    forAll(dpYsource(),ns)
    {
        dpYsource()[ns] = 0.0;
    }
    
    forAll(dpXiCsource(),nXi)
    {
        dpXiCsource()[nXi] = 0.0;
    }
    
    dp_hAsource() = 0.0;

    scalar sumY = 0.0;
    forAll(this->Y(),speciesI)
       sumY+=this->Y()[speciesI];
    forAll(this->Y(),speciesI)
       this->Y()[speciesI]/=sumY;
}


// ************************************************************************* //
