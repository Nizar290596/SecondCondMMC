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

#ifndef EquilibriumApproximation_H
#define EquilibriumApproximation_H

#include "EnvelopeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class EquilibriumApproximation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EquilibriumApproximation
:
    public EnvelopeModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("equilibriumApproximation");

    // Constructors

        //- Construct from dictionary
        EquilibriumApproximation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        EquilibriumApproximation(const EquilibriumApproximation<CloudType>& ea);

        //- Construct and return a clone
        virtual autoPtr<EnvelopeModel<CloudType> > clone() const
        {
            return autoPtr<EnvelopeModel<CloudType> >
            (
                new EquilibriumApproximation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~EquilibriumApproximation();

    bool printOut_;
    label nElements_;
    label nSpeciesSt_;
    label inertStart_;
    List<bool> inertSpecies_;
    List<labelField> elements_;
    wordList elementNames_;
    List<word> speciesSt_;
    List<scalar> fuelRatio_;
    word oxidizer_;
    scalar oxRatio_;
    scalar MW0_;
    scalarField X0_;
    scalarField n0_;
/*
    scalarField n0_;
    scalarField noFuel_;
    scalarField n1_;
    scalar XOx_;
*/
    scalar nOxNet_;
    bool envelopeSwitch_;

    labelField stThermoId_;

    scalar fSt_;
    scalar fSurfOld_;
    scalarField YSt_;// Member Functions

        // Edit
/*
        inline scalarField& n0();
        inline scalarField& noFuel();
        inline scalarField& n1();
        inline scalarField& XOx();
*/
//        inline scalar& nOxNet();
        
        // Access

        //-Read information from input
        void setTableInfo();

        //- initialize parcel
        virtual void initialize
        (
            const scalarField Y0,
            const scalar T0,
	    const label sCell
        );

        //-Gas calculation
        void gasTable
        (
            const scalarField Y0
        );

        //-Vapor calculation
        void solveTable
        (
            const scalarField Y1
        );

        //-Surface mixture fraction
//**//        virtual scalarField fSurface
        virtual scalar fSurface //**//
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar f0,
            const scalarField YSurf
        );
/*
	//- mixture fraction: no envelope
        virtual scalarField fNoEnvelope
        (
            const scalarField YSt,
            const scalarField Y1,
            const scalar fSt,
            const scalarField YX
        );
*/
	//-Surface composition
        virtual scalarField YSurface
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar f0,
            const scalar fSurf
        ) const;

	//-Surface composition
        virtual scalarField YFilm
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar f0,
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

        scalar haEnvelope
        (
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

	scalar TEnvelope
        (
            const scalarField YSurf,
            const scalar haSurf,
            const scalar p,
            const scalar T0
        ) const;

	scalarField YEnvelope
        (
           const scalarField Y0,
           const scalarField Y1,
           const scalar f0,
	   const scalar f1,
           const scalar fX
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EquilibriumApproximation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
