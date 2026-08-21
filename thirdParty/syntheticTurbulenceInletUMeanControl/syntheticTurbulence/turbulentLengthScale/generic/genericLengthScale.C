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
    Foam::genericLengthScale

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

#include "genericLengthScale.H"

namespace Foam
{
	defineTypeNameAndDebug(genericLengthScale, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::genericLengthScale::genericLengthScale
(
	const Time& runTime,
	const dictionary& dict,
	const fvMesh& mesh
)
:
    turbulentLengthScale(),
	LName_(dict.lookup("LName")),
	L_
	(
		mesh.C().size(),
		-1
	)
{

	volScalarField* L = new volScalarField
	(
		IOobject
		(
			LName_,
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	regIOobject::store(L);

	L_ = L->internalField();

}


Foam::genericLengthScale::genericLengthScale
(
	const Time& runTime,
	const dictionary& dict,
	const fvPatchVectorField& pvf
)
:
    turbulentLengthScale(),
	LName_(dict.lookup("LName")),
	L_
	(
		pvf.size(),
		-1
	)
{

	const fvMesh& mesh = pvf.patch().boundaryMesh().mesh();

	volScalarField* L = new volScalarField
	(
		IOobject
		(
			LName_,
			runTime.timeName(),
			runTime,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);
	
	regIOobject::store(L);
	
	label thisPatchID = pvf.patch().index();

	L_ = L->boundaryField()[thisPatchID];

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::genericLengthScale::~genericLengthScale()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField& Foam::genericLengthScale::L() const
{
	return L_;
}


void Foam::genericLengthScale::write(Ostream& os)
{

	os.writeKeyword("L") 	  << nl	<< token::BEGIN_BLOCK << nl;
	os.writeKeyword("type") 	  << typeName << token::END_STATEMENT << nl;
	os.writeKeyword("LName") 	  << LName_ << token::END_STATEMENT << nl;
	os << token::END_BLOCK << endl;	
}

// ************************************************************************* //
