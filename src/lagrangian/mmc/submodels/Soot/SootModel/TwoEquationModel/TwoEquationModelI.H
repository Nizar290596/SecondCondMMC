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

// * * * * * * * * * * * * * * *Member Functions * * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::wordList& Foam::TwoEquationModel<CloudType>::species() const
{
    return this->species_;
}

template<class CloudType>
inline const Foam::scalar& Foam::TwoEquationModel<CloudType>::rhoSoot() const
{
    return this->rhoSoot_;
}

template<class CloudType>
inline const Foam::scalar& Foam::TwoEquationModel<CloudType>::nSootScaling() const
{
    return this->nSootScaling_;
}

template<class CloudType>
inline const Foam::scalar& Foam::TwoEquationModel<CloudType>::surfSoot() const
{
    return this->surfSoot_;
}

template<class CloudType>
inline const Foam::List<Foam::Tuple2<Foam::string, Foam::scalarList> >& Foam::TwoEquationModel<CloudType>::kineticsCoeffs() const
{
    return this->kineticsCoeffs_;
}

template<class CloudType>
inline const Foam::label& Foam::TwoEquationModel<CloudType>::inceptionCarbon() const
{
    return this->inceptionCarbon_;
}

template<class CloudType>
inline const Foam::Switch& Foam::TwoEquationModel<CloudType>::gasSource() const
{
    return this->gasSource_;
}

template<class CloudType>
inline const Foam::List<Foam::Tuple2<Foam::string, Foam::scalar> >& Foam::TwoEquationModel<CloudType>::gasSourceTable() const
{
    return this->gasSourceTable_;
}

template<class CloudType>
inline const Foam::List<Foam::Tuple2<Foam::string, Foam::scalar> >& Foam::TwoEquationModel<CloudType>::YsourceTable() const
{
    return this->YsourceTable_;
}

template<class CloudType>
inline const Foam::scalar& Foam::TwoEquationModel<CloudType>::sootLimit() const
{
    return this->sootLimit_;
}

template<class CloudType>
inline Foam::scalar& Foam::TwoEquationModel<CloudType>::diaSoot()
{
    return this->diaSoot_;
}

template<class CloudType>
inline Foam::scalar& Foam::TwoEquationModel<CloudType>::surfSoot()
{
    return this->surfSoot_;
}

template<class CloudType>
Foam::Switch Foam::TwoEquationModel<CloudType>::turbMixing() const
{
    return this->turbMixing_;
}

template<class CloudType>
inline const Foam::scalarList& Foam::TwoEquationModel<CloudType>::sootContribs() const
{
    return this->sootContribs_;
}

template<class CloudType>
inline Foam::scalarList& Foam::TwoEquationModel<CloudType>::sootContribs()
{
    return this->sootContribs_;
}

// ************************************************************************* //
