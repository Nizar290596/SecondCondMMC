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

    The NoEnvelope model only calls the noEnvelope functions defined in the base class

\*---------------------------------------------------------------------------*/

#include "NoEnvelope.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoEnvelope<CloudType>::NoEnvelope
(
    const dictionary& dict,
    CloudType& owner
)
:
    EnvelopeModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoEnvelope<CloudType>::NoEnvelope
(
    const NoEnvelope<CloudType>& ne
)
:
    EnvelopeModel<CloudType>(ne)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoEnvelope<CloudType>::~NoEnvelope()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
//**//Foam::scalarField Foam::NoEnvelope<CloudType>::fSurface
Foam::scalar Foam::NoEnvelope<CloudType>::fSurface //**//
(
    const scalarField Y0,
    const scalarField Y1,
    const scalar f0,
    const scalarField YSurf
)
{
//**//    scalarField fSurf;
//**//    fSurf.setSize(this->owner().linkFG().size(),-2.0);
    scalar fSurf = -2.0; //**//
    this->setBoundaries(f0,1.0);
    fSurf = this->fNoEnvelope(Y0,Y1,YSurf);

    return fSurf;
}

template<class CloudType>
Foam::scalarField Foam::NoEnvelope<CloudType>::YSurface
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar fSurf
)
{
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YSurf = this->YNoEnvelope(Y0,Y1,fSurf);

    return YSurf;
}

template<class CloudType>
Foam::scalarField Foam::NoEnvelope<CloudType>::YFilm
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar fFilm
) const
{   
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YFilm = this->YNoEnvelope(Y0,Y1,fFilm);
    
    return YFilm;
}

template<class CloudType>
Foam::scalar Foam::NoEnvelope<CloudType>::haSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar ha0,
    const scalar ha1,
    const scalar fSurf
) const
{
    scalar haSurf = -2.0;
    haSurf = this->haNoEnvelope(YSurf,p,TD);

    return haSurf;
}

//-Surface temperature
template<class CloudType>
Foam::scalar Foam::NoEnvelope<CloudType>::TSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar TGas,
    const scalar fSurf,
    const scalar haSurf
) const
{
    scalar TSurf = -2.0;
    TSurf = this->TNoEnvelope(TD);

    return TSurf;
}

// ************************************************************************* //
