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
    Foam::sootPopeCloud

Description
    Virtual abstract base class for templated ItoPopeCloud

SourceFiles
    sootPopeCloud.H

\*---------------------------------------------------------------------------*/

#ifndef sootPopeCloud_H
#define sootPopeCloud_H

#include "typeInfo.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sootPopeCloud Declaration
\*---------------------------------------------------------------------------*/

class sootPopeCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        sootPopeCloud(const sootPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const sootPopeCloud&);


public:

    //- Runtime type information
    TypeName("sootPopeCloud");

    // Constructors

        //- Null constructor
        sootPopeCloud();


    // Member functions


    //- Destructor
    virtual ~sootPopeCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
