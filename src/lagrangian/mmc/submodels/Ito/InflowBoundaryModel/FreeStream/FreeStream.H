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
    Foam::FreeStream

Description
    Inserting new particles across the faces of a all patched of type
    "patch" for a free stream.  Uniform values number density, temperature
    and velocity sourced face-by-face from the boundaryT and boundaryU fields
    of the cloud.

\*---------------------------------------------------------------------------*/

#ifndef FreeStream_H
#define FreeStream_H

#include "InflowBoundaryModel.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class FreeStream Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FreeStream
:
    public InflowBoundaryModel<CloudType>
{
    // Private data

        //- The indices of boundary patches where Pope Particles are introduced
        labelList patches_;

        scalar Cm_;

        //- A List of Lists of Fields specifying carry-over of mass flux from one timestep to the next
        // + Outer List - one inner List for each patch
        // + Inner List - one Field for every species to be introduced
        // + Each field entry corresponding to a face to be injected across
        List<Field<scalar> >  particleFluxAccumulators_;

        List<Field<scalar> >  faceRandThreshold_;

    public:
        //- Runtime type information
        TypeName("FreeStream");

        // Constructors
        //- Construct from dictionary
        FreeStream
        (
            const dictionary& dict,
            CloudType& cloud
        );

        //- Construct copy
        FreeStream(const FreeStream<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InflowBoundaryModel<CloudType> > clone() const
        {
            return autoPtr<InflowBoundaryModel<CloudType> >
            (
                new FreeStream<CloudType>(*this)
            );
        }

        // Destructor
        virtual ~FreeStream();

        // Member Functions

        virtual void inflow();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FreeStream.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


