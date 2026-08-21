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

Class
    Foam::interpolated

Description
    This class is used for the evolution of the referece variable from 
    LES fields.

SourceFiles
    interpolated.C

\*---------------------------------------------------------------------------*/

#ifndef interpolated_H
#define interpolated_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"

#include "evolModel.H"
#include "interpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolated Declaration
\*---------------------------------------------------------------------------*/

class interpolated
:
    public evolModel
{
    protected:

        // Protected data

        word refName_;

//        volScalarField& XiRField_;

//        autoPtr<interpolation<scalar>> XiInterp_;

    public:

        //- Runtime type information
        TypeName("interpolated");

        // Constructors

            //- Construct from dictionary and mesh
            interpolated
            (
                const word& refName,
                const fvMesh& mesh
            );

            //- Construct copy
            interpolated(const interpolated& cm);


        //- Destructor
        virtual ~interpolated()
        {};

        // Member functions

            //- virtual functions
                virtual scalar compute
                (
                    const particle& p,
                    scalar dt, 
                    const scalar& x
                ) const;

                virtual scalar initialise
                (
                    label patch,
                    label patchFace,
                    particle& p
                ) const;

                virtual scalar initialise
                (
                    label celli,
                    particle& p
                ) const;

            //- Construct and return a clone to evolModel
                virtual autoPtr<evolModel> clone() const
                {
                    return autoPtr<evolModel>(new interpolated(*this));
                }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
