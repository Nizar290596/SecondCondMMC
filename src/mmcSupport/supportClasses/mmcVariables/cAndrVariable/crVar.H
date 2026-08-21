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
    Foam::crVar

Description
    This classes is aimed to help the creation of variables that work as
    reference and coupling variables in MMC. It may olso be used for 
    variables required in postprocessing.

SourceFiles
    crVar.C

\*---------------------------------------------------------------------------*/

#ifndef crVar_H
#define crVar_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"

#include "couplingVar.H"
#include "referenceVar.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class crVar Declaration
\*---------------------------------------------------------------------------*/

class crVar
:
    //public couplingVar,
    public referenceVar,
    public couplingVar
{
    public:
    
        //- Runtime type information
        TypeName("couplingAndreferenceVar");
        
        // Constructors
                       
            //- Construct from dictionary and mesh
            crVar
            (
                const dictionary& entry,
                const fvMesh& mesh
            );
            
            //- Construct copy
            crVar(const crVar& cm);//!!


        //- Destructor
        virtual ~crVar()
        {};
            
        // Member functions
        virtual autoPtr<basicMMCVar> clone() const
	{
	    return autoPtr<basicMMCVar>(new crVar(*this));
	}   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
