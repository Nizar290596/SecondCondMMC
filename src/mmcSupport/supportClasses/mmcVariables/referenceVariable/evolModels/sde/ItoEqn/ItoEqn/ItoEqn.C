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

#include "ItoEqn.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class T>
Foam::ItoEqn<T>::ItoEqn()
:
    A_(),
    D_(),
    rndGen_(Pstream::myProcNo()+1)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::ItoEqn<T>::computeDriftCoeff
(
    const particle& p, 
    const T& x
) const
{
    notImplemented
    (
        "void Foam::ItoEqn<T>::computeDriftCoeff"
        "("
            "const particle& p,"
            "const T& x"
        ") const "
    );
}


template<class T>
void Foam::ItoEqn<T>::computeDiffusionCoeff(const particle& p) const
{
    notImplemented
    (
        "void Foam::ItoEqn<T>::computeDiffusionCoeff"
        "("
            "const particle& p,"
        ") const "
    );
}


template<class T>
T Foam::ItoEqn<T>::integrate
(
    const particle& p,
    scalar dt, 
    const T& x
) const
{
    notImplemented
    (
        "T Foam::ItoEqn<T>::integrate"
        "("
            "const particle& p,"
            "scalar dt,"
            "const T& x"
        ") const "
    );

    return 0;
}
// * * * * * * * * * * * * * * Specialized Members * * * * * * * * * * * * * //

#include "scalarItoEqn.C"

// * * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * //

#include "ItoEqnNew.C"

// ************************************************************************* //
