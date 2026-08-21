/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This file is not part of OpenFOAM but an original routine developed 
    in the OpenFOAM techonology.

    You can redistribute it and/or modify it under the terms of the GNU General
    Public License as published by the Free Software Foundation; either version
    3 of the License, or (at your option) any later version. See GNU General
    Public License at <http://www.gnu.org/licenses/gpl.html>

    OpenFOAM is a trademark of OpenCFD.

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
    Foam::turbulentLengthScale

Authors
    Andrea Montorfano, Federico Piscaglia
    Dipartimento di Energia, Politecnico di Milano
    via Lambruschini 4
    I-20156 Milano (MI)
    ITALY

Contact
    andrea.montorfano@polimi.it , ph. +39 02 2399 3909
    federico.piscaglia@polimi.it, ph. +39 02 2399 8620    

\*---------------------------------------------------------------------------*/

#include "turbulentLengthScale.H"

#include "uniformLengthScale.H"
#include "genericLengthScale.H"
#include "kwLengthScale.H"

namespace Foam
{
    defineTypeNameAndDebug(turbulentLengthScale, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentLengthScale::turbulentLengthScale()
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::turbulentLengthScale*
Foam::turbulentLengthScale::New
(
    const dictionary& dict,
    const Time& runTime,
    const Type& support,
  const word& turbulentLengthScaleName
)
{

    if (turbulentLengthScaleName == uniformLengthScale::typeName)
    {
        return     new uniformLengthScale(dict, support);
    }
    else if (turbulentLengthScaleName == kwLengthScale::typeName)
    {
        return new kwLengthScale(runTime, dict, support);
    }
    else if (turbulentLengthScaleName == genericLengthScale::typeName)
    {
        return new genericLengthScale(runTime, dict, support);
    }
    else
    {
        FatalErrorIn("Foam::turbulentLengthScale::New")
            << "Unknown length scale model " << turbulentLengthScaleName
            << exit(FatalError);
    }

    return NULL;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentLengthScale::~turbulentLengthScale()
{}

template  Foam::turbulentLengthScale* Foam::turbulentLengthScale::New(const dictionary&, const Time&, const fvPatchVectorField&, const word&);
template  Foam::turbulentLengthScale* Foam::turbulentLengthScale::New(const dictionary&, const Time&, const fvMesh&, const word&);
// ************************************************************************* //
