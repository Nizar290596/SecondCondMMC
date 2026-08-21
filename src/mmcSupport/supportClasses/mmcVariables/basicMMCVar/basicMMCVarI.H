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

inline const Foam::word& Foam::basicMMCVar::name() const
{
    return name_;
}


inline Foam::volScalarField& Foam::basicMMCVar::field()
{
    return field_();
}


inline const Foam::volScalarField& Foam::basicMMCVar::field() const
{
    return field_();
}


inline Foam::word Foam::basicMMCVar::cName() const
{
    return name_;
}


inline Foam::word Foam::basicMMCVar::rName() const
{
    return name_;
}


inline Foam::word Foam::basicMMCVar::refType() const
{
    return word::null;
}


inline Foam::Switch Foam::basicMMCVar::passive() const
{
    return true;
}


inline Foam::label Foam::basicMMCVar::mixStep() const
{
    return -1;// If mixingin step = -1 don't mix
}

inline const Foam::evolModel& Foam::basicMMCVar::evolveMethod() const
{
    return evolModel::null();
}
// ************************************************************************* //
