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


template<class CloudType, class ReactionThermo>
inline const ReactionThermo& Foam::ThermoPopeCloud<CloudType,ReactionThermo>::thermo() const
{
    return *thermo_;
}


template<class CloudType, class ReactionThermo>
inline ReactionThermo& Foam::ThermoPopeCloud<CloudType,ReactionThermo>::thermo()
{
    return *thermo_;
}


template<class CloudType, class ReactionThermo>
inline const Foam::SLGThermo& Foam::ThermoPopeCloud<CloudType,ReactionThermo>::slgThermo() const
{
    return SLGThermo_();
}


template<class CloudType, class ReactionThermo>
inline const Foam::volScalarField& Foam::ThermoPopeCloud<CloudType,ReactionThermo>::T() const
{
    return T_;
}


template<class CloudType, class ReactionThermo>
inline const Foam::volScalarField& Foam::ThermoPopeCloud<CloudType,ReactionThermo>::p() const
{
    return p_;
}


template<class CloudType, class ReactionThermo>
inline const Foam::CompositionModel<Foam::ThermoPopeCloud<CloudType,ReactionThermo> >&
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::composition() const
{
    return compositionModel_();
}


template<class CloudType, class ReactionThermo>
inline const Foam::ThermoPhysicalCouplingModel<Foam::ThermoPopeCloud<CloudType,ReactionThermo> >&
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::coupling() const
{
    return thermoPhysicalCouplingModel_();
}


template<class CloudType, class ReactionThermo>
inline Foam::CompositionModel<Foam::ThermoPopeCloud<CloudType,ReactionThermo> >&
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::composition() 
{
    return compositionModel_();
}


// - Return access to density coupling model
template<class CloudType, class ReactionThermo>
inline Foam::ThermoPhysicalCouplingModel<Foam::ThermoPopeCloud<CloudType,ReactionThermo> >&
Foam::ThermoPopeCloud<CloudType,ReactionThermo>::coupling()
{
    return thermoPhysicalCouplingModel_();
}


template<class CloudType, class ReactionThermo>
inline Foam::scalar Foam::ThermoPopeCloud<CloudType,ReactionThermo>::Tmax() const
{
    scalar T = -GREAT;
    scalar n = 0;
    forAllConstIters(*this, iter)
    {
        const particleType& p = iter();
        T = max(T, p.T());
        n++;
    }

    reduce(T, maxOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType, class ReactionThermo>
inline Foam::scalar Foam::ThermoPopeCloud<CloudType,ReactionThermo>::Tmin() const
{
    scalar T = GREAT;
    scalar n = 0;
    forAllConstIters(*this, iter)
    {
        const particleType& p = iter();
        T = min(T, p.T());
        n++;
    }

    reduce(T, minOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}

// ************************************************************************* //
