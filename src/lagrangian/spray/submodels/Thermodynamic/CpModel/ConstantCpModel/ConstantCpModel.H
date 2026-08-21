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
    Foam::ConstantCpModel

Description
    Inserting new particles across the faces of a all patched of type
    "patch" for a free stream.  Uniform values number density, temperature
    and velocity sourced face-by-face from the boundaryT and boundaryU fields
    of the cloud.

\*---------------------------------------------------------------------------*/

#ifndef ConstantCpModel_H
#define ConstantCpModel_H

#include "CpModel.H"
#include "List.H"
#include "scalar.H"
#include "scalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class ConstantCpModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConstantCpModel
:
    public CpModel<CloudType>
{
    // Private data
        // Constant Specific Heat of Dispersed Phase
        scalar Cp_;

    public:
        //- Runtime type information
        TypeName("ConstantCp");

    // Constructors

        //- Construct from dictionary
        ConstantCpModel(const dictionary&, CloudType& owner);

        //- Construct copy
        ConstantCpModel(const ConstantCpModel<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<CpModel<CloudType>> clone() const
        {
            return autoPtr<CpModel<CloudType>>
            (
                new ConstantCpModel<CloudType>(*this)
            );
        }

        // Destructor
        virtual ~ConstantCpModel();


        // Member Functions

        // Specific heat at constant pressure of Dispersed Phase
        virtual scalar Cp(const scalar) const;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ConstantCpModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


