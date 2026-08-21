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
    Foam::uniformLengthScale


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

#include "uniformLengthScale.H"

namespace Foam
{
	defineTypeNameAndDebug(uniformLengthScale, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformLengthScale::uniformLengthScale
(
	const dictionary& dict,
	const fvMesh& mesh
)
:
    turbulentLengthScale(),
	L_
	(
		mesh.C().size(),
		readScalar(dict.lookup("value"))
	)
{}


Foam::uniformLengthScale::uniformLengthScale
(
	const dictionary& dict,
	const fvPatchVectorField& pvf
)
:
    turbulentLengthScale(),
	L_
	(
		pvf.size(),
		readScalar(dict.lookup("value"))
	)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformLengthScale::~uniformLengthScale()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField& Foam::uniformLengthScale::L() const
{
	return L_;
}


void Foam::uniformLengthScale::write(Ostream& os)
{
	os.writeKeyword("L") 	  << nl	<< token::BEGIN_BLOCK << nl;
	os.writeKeyword("type") 	  << typeName << token::END_STATEMENT << nl;
	os.writeKeyword("value") 	  << average(L_) << token::END_STATEMENT << nl;
	os << token::END_BLOCK << endl;	
}

// ************************************************************************* //
