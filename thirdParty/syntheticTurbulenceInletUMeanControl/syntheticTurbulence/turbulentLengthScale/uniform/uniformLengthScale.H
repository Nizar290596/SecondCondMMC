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

Description
    Support class for synthetic turbulence generation. Specifies an uniform
    turbulent integral scale over the support entity (either a boundary patch
    or a volume).
    
    Example:
    
    \verbatim
    L
    {
        type    uniformLengthScale;
        value   0.05;
    }
    \endverbatim
    
Authors
	Andrea Montorfano, Federico Piscaglia
	Dipartimento di Energia, Politecnico di Milano
	via Lambruschini 4
	I-20156 Milano (MI)
	ITALY

Contact
	andrea.montorfano@polimi.it , ph. +39 02 2399 3909
	federico.piscaglia@polimi.it, ph. +39 02 2399 8620	

SourceFiles
    uniformLengthScale.C

\*---------------------------------------------------------------------------*/

#ifndef uniformLengthScale_H
#define uniformLengthScale_H

#include "turbulentLengthScale.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class uniformLengthScale Declaration
\*---------------------------------------------------------------------------*/

class uniformLengthScale
:
    public turbulentLengthScale
{
    // Private data

        //- Description of data_
        scalarField L_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformLengthScale(const uniformLengthScale&);

        //- Disallow default bitwise assignment
        void operator=(const uniformLengthScale&);


public:

    //- Runtime type information
    TypeName("uniformLengthScale");


    // Constructors

        //- Construct from fvMesh
        uniformLengthScale
		(
			const dictionary& dict, 
			const fvMesh& mesh
		);

        //- Construct from fvPatchVectorField
        uniformLengthScale
		(
			const dictionary& dict, 
			const fvPatchVectorField& pvf
		);

    //- Destructor
    virtual ~uniformLengthScale();


    // Member Functions

        //- return length
		virtual const scalarField& L() const;

        //- write output
		virtual void write(Ostream& os);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
