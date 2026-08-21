/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
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

#ifndef NoMixing_H
#define NoMixing_H

#include "CloudMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoMixing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoMixing
:
    public CloudMixingModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from dictionary
        NoMixing(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        NoMixing(const NoMixing<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new NoMixing<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoMixing();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoMixing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
