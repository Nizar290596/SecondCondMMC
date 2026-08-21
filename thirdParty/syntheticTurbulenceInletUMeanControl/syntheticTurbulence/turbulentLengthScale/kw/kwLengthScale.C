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
    Foam::kwLengthScale


Authors
	Andrea Montorfano, Federico Piscaglia
	Dipartimento di Energia, Politecnico di Milano
	via Lambruschini 4
	I-20156 Milano (MI)
	ITALY

Contact
	andrea.montorfano@polimi.it , ph. +39 02 2399 3909
	federico.piscaglia@polimi.it, ph. +39 02 2399 8620	

\* ---------------------------------------------------------------------------*/

#include "kwLengthScale.H"
#include "objectRegistry.H"
#include "fvPatch.H"

namespace Foam
{
	defineTypeNameAndDebug(kwLengthScale, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kwLengthScale::kwLengthScale
(
	const Time& runTime,
	const dictionary& dict,
	const fvMesh& mesh
)
:
    turbulentLengthScale(),
	L_
	(
		mesh.C().size(),
		-1
	),
	kName_(dict.lookupOrDefault<word>("kName","k")),
	omegaName_(dict.lookupOrDefault<word>("omegaName","omega")),
	fL_(1.0)
{

	if (dict.found("fL"))
	{
		fL_ = readScalar(dict.lookup("fL"));
	}

	const volScalarField& k = runTime.lookupObject<volScalarField>(kName_);
	
	const volScalarField& omega = runTime.lookupObject<volScalarField>(omegaName_);

	volScalarField* L = new volScalarField
	(
		IOobject
		(
			"L",
			runTime.timeName(),
			runTime,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fL_*sqrt(k)/omega
	);
	
	regIOobject::store(L);


	L_ = L->internalField();

}


Foam::kwLengthScale::kwLengthScale
(
	const Time& runTime,
	const dictionary& dict,
	const fvPatchVectorField& pvf
)
:
    turbulentLengthScale(),
	L_
	(
		pvf.size(),
		-1
	),
	kName_(dict.lookupOrDefault<word>("kName","k")),
	omegaName_(dict.lookupOrDefault<word>("omegaName","omega")),
	fL_(1.0)
{

	if (dict.found("fL"))
	{
		fL_ = readScalar(dict.lookup("fL"));
	}

	Info << "kwLengthScale" << endl;
	//const fvMesh& mesh = pvf.patch().boundaryMesh().mesh();

	const volScalarField& k = runTime.lookupObject<volScalarField>(kName_);
	const volScalarField& omega = runTime.lookupObject<volScalarField>(omegaName_);

	label thisPatchID = pvf.patch().index();

	const fvPatchScalarField& kBnd = k.boundaryField()[thisPatchID];
	const fvPatchScalarField& omegaBnd = omega.boundaryField()[thisPatchID];	

	L_ = fL_*sqrt(kBnd)/omegaBnd;

	//- store volScalarField for later use or debug 
	/*volScalarField* L = new volScalarField
	(
		IOobject
		(
			"L",
			runTime.timeName(),
			runTime,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("L",dimLength,pTraits<scalar>::zero),
		k.boundaryField().types()
	);

	L->boundaryField()[thisPatchID] = L_;

	regIOobject::store(L);*/
	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::kwLengthScale::~kwLengthScale()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField& Foam::kwLengthScale::L() const
{
	return L_;
}


void Foam::kwLengthScale::write(Ostream& os)
{
	os.writeKeyword("L") 	  		<< nl	<< token::BEGIN_BLOCK << nl;
	os.writeKeyword("type") 	  	<< typeName << token::END_STATEMENT << nl;
	os.writeKeyword("kName") 	  	<< kName_ << token::END_STATEMENT << nl;
	os.writeKeyword("omegaName") 	<< omegaName_ << token::END_STATEMENT << nl;	
	os.writeKeyword("fL") 	  		<< fL_ << token::END_STATEMENT << nl;	
	os << token::END_BLOCK << endl;	
}

// ************************************************************************* //
