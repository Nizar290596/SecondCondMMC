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
    Foam::NoInflow

Description
    Not inserting any particles

\*---------------------------------------------------------------------------*/

#ifndef NoInflow_H
#define NoInflow_H

#include "InflowBoundaryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class NoInflow Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoInflow
:
    public InflowBoundaryModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("NoInflow");


    // Constructors

        //- Construct from dictionary
        NoInflow
        (
            const dictionary& dict,
            CloudType& cloud
        );
        //- Construct copy
        NoInflow(const NoInflow<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InflowBoundaryModel<CloudType> > clone() const
        {
            return autoPtr<InflowBoundaryModel<CloudType> >
            (
                new InflowBoundaryModel<CloudType>(*this)
            );
        }

    // Destructor
    virtual ~NoInflow();


    // Member Functions

    virtual void inflow();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoInflow.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
