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

#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline bool Foam::SootPopeCloud<CloudType>::radiationFlag() const
{
    return radiation_;
}


template<class CloudType>
inline const Foam::Switch& Foam::SootPopeCloud<CloudType>::sootingFlame() const
{
    return sootingFlame_;
}


template<class CloudType>
inline const Foam::scalar& Foam::SootPopeCloud<CloudType>::sootSolidDensity() const
{
    return sootSolidDensity_;
}


template<class CloudType>
inline const Foam::List<Foam::scalar>& Foam::SootPopeCloud<CloudType>::sootMinMolWt() const
{
    return sootMinMolWt_;
}


template<class CloudType>
inline const Foam::label& Foam::SootPopeCloud<CloudType>::numfVsoot() const
{
    return numfVsoot_;
}


template<class CloudType>
inline const Foam::scalar& Foam::SootPopeCloud<CloudType>::sootMinKn() const
{
    return sootMinKn_;
}


template<class CloudType>
inline const Foam::scalar& Foam::SootPopeCloud<CloudType>::IntermittencyfV() const
{
    return IntermittencyfV_;
}


template<class CloudType>
inline const Foam::SootModel<Foam::SootPopeCloud<CloudType>>&
Foam::SootPopeCloud<CloudType>::soot() const 
{
    return sootModel_();
}


template<class CloudType>
inline const Foam::CloudRadiationModel<Foam::SootPopeCloud<CloudType> >&
Foam::SootPopeCloud<CloudType>::radiation() const
{
    return radiationModel_();
}



template<class CloudType>
inline Foam::SootModel<Foam::SootPopeCloud<CloudType>>&
Foam::SootPopeCloud<CloudType>::soot()
{
    return sootModel_();
}


template<class CloudType>
inline Foam::CloudRadiationModel<Foam::SootPopeCloud<CloudType> >&
Foam::SootPopeCloud<CloudType>::radiation() 
{
    return radiationModel_();
}


template<class CloudType>
inline const Foam::List<Foam::label>& 
Foam::SootPopeCloud<CloudType>::indexO2N2() const
{
    return indexO2N2_;
}


// ************************************************************************* //
