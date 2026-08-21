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
    Foam::ThinRadiationModel

\*---------------------------------------------------------------------------*/

#ifndef ThinRadiationModel_H
#define ThinRadiationModel_H

#include "FPRadiationModel.H"
#include "List.H"
#include "scalar.H"
#include "scalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class ThinRadiationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThinRadiationModel
:
    public FPRadiationModel<CloudType>
{
    // Private data
        // Reference/Ambient Temperature
        scalar TRef_;


    public:
        //- Runtime type information
        TypeName("OpticallyThinRadiation");

    // Constructors

        //- Construct from dictionary
        ThinRadiationModel(const dictionary&, CloudType& owner);

        //- Construct copy
        ThinRadiationModel(const ThinRadiationModel<CloudType>& trm);

        //- Construct and return a clone
        virtual autoPtr<FPRadiationModel<CloudType>> clone() const
        {
            return autoPtr<FPRadiationModel<CloudType>>
            (
                new ThinRadiationModel<CloudType>(*this)
            );
        }

        // Destructor
        virtual ~ThinRadiationModel();


        // Member Functions

        // Radiative Flux
        virtual scalar qR
        (
        	const scalar TD,   // Droplet Temperature
            const scalar mFlux // Evaporative mass flux
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThinRadiationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
