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

Description
    Abstract base class for turbulence integral length scale.
  


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
    turbulentLengthScale.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentLengthScale_H
#define turbulentLengthScale_H

#include "fvCFD.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class turbulentLengthScale Declaration
\*---------------------------------------------------------------------------*/


class turbulentLengthScale
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulentLengthScale(const turbulentLengthScale&);

        //- Disallow default bitwise assignment
        void operator=(const turbulentLengthScale&);

public:

    //- Runtime type information
    TypeName("turbulentLengthScale");

    // Constructors

        //- construct null
        turbulentLengthScale();

    // Selectors

        //- Return a pointer to the selected TLS model,
        // to use either as volField or patchField
        template<class Type>
        static turbulentLengthScale* New (
                const dictionary& dict,
                const Time& runTime,
                const Type& support,
          const word& turbulentLengthScaleName = turbulentLengthScale::typeName
        );


    //- Destructor
    virtual ~turbulentLengthScale();

    // Member Functions

        //- return length
        virtual const scalarField& L() const=0;
        

        //- Write output
        virtual void write(Ostream& os)=0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



// ************************************************************************* //
