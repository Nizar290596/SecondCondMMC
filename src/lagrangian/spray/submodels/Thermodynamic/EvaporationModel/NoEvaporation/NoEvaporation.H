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

#ifndef NoEvaporation_H
#define NoEvaporation_H

#include "EvaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoEvaporation
:
    public EvaporationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from dictionary
        NoEvaporation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoEvaporation(const NoEvaporation<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<EvaporationModel<CloudType> > clone() const
        {
            return autoPtr<EvaporationModel<CloudType> >
            (
                new NoEvaporation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoEvaporation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
