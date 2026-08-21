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
Class
    Foam::VariableHardSphere

Description
    Variable Hard Sphere BinaryCollision Model

\*---------------------------------------------------------------------------*/

#ifndef NoEnvelope_H
#define NoEnvelope_H

#include "EnvelopeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoEnvelope Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoEnvelope
:
    public EnvelopeModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from dictionary
        NoEnvelope(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoEnvelope(const NoEnvelope<CloudType>& ne);

        //- Construct and return a clone
        virtual autoPtr<EnvelopeModel<CloudType> > clone() const
        {
            return autoPtr<EnvelopeModel<CloudType> >
            (
                new NoEnvelope<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoEnvelope();

    // Member Functions

        // Access


        //-Surface mixture fraction
//**//        virtual scalarField fSurface
        virtual scalar fSurface //**//
        (
            const scalarField YFuelVap,
            const scalarField Y1,
            const scalar f0,
            const scalarField YSurf
        );

	//-Surface composition
        virtual scalarField YSurface
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar fSurf
        );

	//-Surface composition
        virtual scalarField YFilm
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar fFim
        ) const;

	//-Surface total enthalpy
        virtual scalar haSurface
        (
            const scalarField YSurf,
            const scalar p,
            const scalar Td,
	        const scalar ha0,
	        const scalar ha1,
	        const scalar fSurf
        ) const;

        //-Surface temperature
        virtual scalar TSurface
        (
            const scalarField YSurf,
            const scalar p,
            const scalar Td,
            const scalar TGas,
            const scalar fSurf,
            const scalar haSurf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoEnvelope.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
