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
    Foam::popeCloudSolution

Description
    Stores all relevant solution info for cloud

SourceFiles
    popeCloudSolutionI.H
    popeCloudSolution.C

\*---------------------------------------------------------------------------*/

#ifndef popeCloudSolution_H
#define popeCloudSolution_H

#include "fvMesh.H"
#include "Switch.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class popeCloudSolution Declaration
\*---------------------------------------------------------------------------*/

class popeCloudSolution
{
    // Private Data

        //- Dictionary used during construction
        dictionary dict_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const popeCloudSolution&);


public:

    // Constructors

        //- Construct from dictionary
        popeCloudSolution(const dictionary& dict);

        //- Construct copy
        popeCloudSolution(const popeCloudSolution& cs);


    //- Destructor
    virtual ~popeCloudSolution();


    // Member functions

        //- Read properties from dictionary
        void read();


        // Access

            //- Return const access to the dictionary
            inline const dictionary& dict() const;

            //- Interpolation schemes dictionary
            inline const dictionary& interpolationSchemes() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "popeCloudSolutionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
