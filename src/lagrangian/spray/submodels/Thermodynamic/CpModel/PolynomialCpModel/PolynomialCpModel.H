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

Description
    An more accurate specific heat capacity (Cp) model using five 
    constants. 

\*---------------------------------------------------------------------------*/

#ifndef PolynomialCpModel_H
#define PolynomialCpModel_H

#include "FixedList.H"
#include "CpModel.H"
#include "scalar.H"
#include "scalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class PolynomialCp Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PolynomialCpModel
:
    public CpModel<CloudType>
{

public:

    // Public data

        static const int nCoeffs_ = 5;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private data
        coeffArray CpCoeffs_;

public:
    //- Runtime type information
    TypeName("PolynomialCp");

 
    // Constructors

        //- Construct from dictionary
        PolynomialCpModel(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        PolynomialCpModel(const PolynomialCpModel<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<CpModel<CloudType>> clone() const
        {
            return autoPtr<CpModel<CloudType>>
            (
                new PolynomialCpModel<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~PolynomialCpModel();


    // Member Functions

        // Specific heat at constant pressure of Dispersed Phase
        virtual scalar Cp(const scalar T) const;

        //Access

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& CpCoeffs() const;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "PolynomialCpModelI.H"

#ifdef NoRepository
#   include "PolynomialCpModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


