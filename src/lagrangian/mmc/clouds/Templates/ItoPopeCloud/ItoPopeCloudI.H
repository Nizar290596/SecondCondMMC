/*---------------------------------------------------------------------------*\
 ##   ####  ######     | 
 ##  ##     ##         | Copyright: ICE Stroemungsfoschungs GmbH
 ##  ##     ####       |
 ##  ##     ##         | http://www.ice-sf.at
 ##   ####  ######     |
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ItoPopeCloud

Description
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline void Foam::ItoPopeCloud<CloudType>::storeDeltaT()
{
    cachedDeltaT_ = mesh().time().deltaT().value();
}


template<class CloudType>
inline Foam::scalar Foam::ItoPopeCloud<CloudType>::cachedDeltaT() const
{
    return cachedDeltaT_;
}


template<class CloudType>
inline const Foam::InflowBoundaryModel<Foam::ItoPopeCloud<CloudType> >& Foam::ItoPopeCloud<CloudType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class CloudType>  
inline Foam::InflowBoundaryModel<Foam::ItoPopeCloud<CloudType> > & Foam::ItoPopeCloud<CloudType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::ItoPopeCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline Foam::scalar Foam::ItoPopeCloud<CloudType>::deltaM() const
{
    return deltaM_;
}


template<class CloudType>
inline const Foam::Time & Foam::ItoPopeCloud<CloudType>::runTime() const
{
    return runTime_;
}

 
template<class CloudType>
inline Foam::StochasticLib1& Foam::ItoPopeCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ItoPopeCloud<CloudType>::delta() const
{
   return turbModelPtr_->delta();
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ItoPopeCloud<CloudType>::muEff() const
{
    return turbModelPtr_->mu() + turbModelPtr_->mut();
}  
  
template<class CloudType>
inline const Foam::volScalarField& Foam::ItoPopeCloud<CloudType>::kfv() const
{
   return turbModelPtr_->k();
}  

template<class CloudType>
inline const Foam::volVectorField& Foam::ItoPopeCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ItoPopeCloud<CloudType>::DEff() const
{
    return DEff_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::ItoPopeCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::ItoPopeCloud<CloudType>::gradRho() const
{
    return gradRho_;
}


template<class CloudType>
inline const Foam::IOdictionary& Foam::ItoPopeCloud<CloudType>::cloudProperties() const
{
    return cloudProperties_;
}


template<class CloudType>
inline const Foam::dictionary& Foam::ItoPopeCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline Foam::popeCloudSolution& Foam::ItoPopeCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const Foam::popeCloudSolution& Foam::ItoPopeCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::eulerianStatistics& Foam::ItoPopeCloud<CloudType>::eulerianStats()
{
    return eulerianStats_;
}


template<class CloudType>
inline const Foam::eulerianStatistics& Foam::ItoPopeCloud<CloudType>::eulerianStats() const
{
    return eulerianStats_;
}


template<class CloudType>
inline const Foam::dictionary& Foam::ItoPopeCloud<CloudType>::eulerianStatsDict() const
{
    return eulerianStatsDict_;
}


template<class CloudType>
inline void Foam::ItoPopeCloud<CloudType>::clear()
{
    return IDLList<typename CloudType::particleType>::clear();
}


template<class CloudType>
inline Foam::scalar& Foam::ItoPopeCloud<CloudType>::fvMinusPpMass()
{
    return fvMinusPpMass_;
}


template<class CloudType>
inline const Foam::scalarField& Foam::ItoPopeCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline const Foam::Switch& Foam::ItoPopeCloud<CloudType>::writeMassConsistency() const
{
    return writeMassConsistency_;
}
// ************************************************************************* //
