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
    Foam::mixingPopeCloud

Description
    Virtual abstract base class for templated ItoPopeCloud

SourceFiles
    mixingPopeCloud.H

\*---------------------------------------------------------------------------*/

#ifndef mixingPopeCloud_H
#define mixingPopeCloud_H

#include "typeInfo.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mixingPopeCloud Declaration
\*---------------------------------------------------------------------------*/

class mixingPopeCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixingPopeCloud(const mixingPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const mixingPopeCloud&);


public:

    //- Runtime type information
    TypeName("mixingPopeCloud");

    // Constructors

        //- Null constructor
        mixingPopeCloud();


    // Member functions

    //- Destructor
    virtual ~mixingPopeCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
