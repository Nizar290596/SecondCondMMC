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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::varSet& Foam::varSet::null()
{
    return *reinterpret_cast< varSet* >(0);
}


inline const Foam::varTable& 
Foam::varSet::varNames() const
{
    return this->varNames_;
}


inline Foam::PtrList<Foam::basicMMCVar>&
Foam::varSet::Vars()
{
    return this->Vars_;
}


inline const Foam::PtrList<Foam::basicMMCVar>&
Foam::varSet::Vars() const
{
    return this->Vars_;
}


inline Foam::basicMMCVar& 
Foam::varSet::Vars(const label i)
{
    return Vars_[i];
}


inline const Foam::basicMMCVar& 
Foam::varSet::Vars(const label i) const
{
    return Vars_[i];
}


inline Foam::basicMMCVar& 
Foam::varSet::Vars(const word& varName)
{
    return Vars_[varNames_[varName]];
}


inline const Foam::basicMMCVar& 
Foam::varSet::Vars(const word& varName) const
{
    return Vars_[varNames_[varName]];
}


inline bool Foam::varSet::contains(const word& varName) const
{
    return varNames_.found(varName);
}

// ************************************************************************* //
