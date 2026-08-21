/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "PolynomialCpModel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PolynomialCpModel<CloudType>::PolynomialCpModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    CpModel<CloudType>(dict,owner,typeName),
    CpCoeffs_(this->coeffDict().lookup("CpDCoeffs"))
{}

template<class CloudType>
Foam::PolynomialCpModel<CloudType>::PolynomialCpModel(const PolynomialCpModel<CloudType>& cpm)
:
    CpModel<CloudType>(cpm),
    CpCoeffs_(this->coeffDict().lookup("CpDCoeffs"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PolynomialCpModel<CloudType>::~PolynomialCpModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
Foam::scalar Foam::PolynomialCpModel<CloudType>::Cp
(
    const scalar T
) const    
{

    const coeffArray& a = CpCoeffs_;
    return
        ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
}

// ************************************************************************* //
