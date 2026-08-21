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
    //- number of elements considered in the table
    label nElements_;

    //- number of species at stoichiometry
    label nSpeciesSt_;
//    label inertStart_;
//    List<bool> inertSpecies_;

    //- species at stoichiometry
    List<bool> stSpecies_;

    //- table of elemnt numbers for each species
    List<labelField> elements_;

    //- names of the considered elements
    wordList elementNames_;

    //- species names of products
    List<word> speciesSt_;

    //- number of elements in products (each elements except the oxidizer element can only be part of one product)
    List<scalar> fuelRatio_;

    //- species name of the oxidizer
    word oxidizer_;

    //- number of elements in oxidizer
    scalar oxRatio_;

    //- Molecular weight of the gas composition
    scalar MW0_;

    //- mole fractions of the gas composition
    scalarField X0_;

    //- number of each element in the gas
    scalarField n0_;

    //- mols of excess oxidizer in gas
    scalar nOxNet_;

    //- switch for envelope state
    bool envelopeSwitch_;

    //- thermo IDs of species at stoichiometry
    labelField stThermoId_;

    //- stoichimetric mixture fraction
    scalar fSt_;

    //- stoichiometry at the surface at the previous iteration
    scalar fSurfOld_;

    //- Fuel (vapor) composition at the previous iteration (only necessary for diff evap)
    scalarField Y1LiqOld_;

    //- Composition at stoichiometry
    scalarField YSt_;// Member Functions

        // Edit
        
        // Access

        //- Read information from input
        void setTableInfo();

        //- initialize parcel
        virtual void initialize
        (
            const scalarField Y0,
            const scalar T0,
	    const label sCell
        );

        //- Gas calculation
        void gasTable
        (
            const scalarField Y0
        );

        //- Vapor calculation
        void solveTable
        (
            const scalarField Y1,
            const scalar f0
        );

        //- Surface mixture fraction
//**//        virtual scalarField fSurface
        virtual scalar fSurface //**//
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar f0,
            const scalarField YSurf
        );

	//- Surface composition
        virtual scalarField YSurface
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar fSurf
        );

	//- Film composition
        virtual scalarField YFilm
        (
            const scalarField Y0,
            const scalarField Y1Liq,
            const scalar fFim
        ) const;

	//- Surface total enthalpy
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

        scalar fEnvelope
        (
            const scalarField Y0,
            const scalarField Y1,
            const scalar f0,
            const scalar f1,
            const scalarField YX
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
