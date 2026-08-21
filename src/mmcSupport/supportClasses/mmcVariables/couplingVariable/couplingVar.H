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
    Foam::couplingVar

Description
    This class is aimed to help the creation of coupling variables in MMC. 
    It may olso be used for variables required in postprocessing.

SourceFiles
    couplingVar.C

\*---------------------------------------------------------------------------*/

#ifndef couplingVar_H
#define couplingVar_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"

#include "basicMMCVar.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class couplingVar Declaration
\*---------------------------------------------------------------------------*/

class couplingVar
:
    public virtual basicMMCVar
{
    protected:
            
        // Protected data
        
            //- coupling name (name in particles)
            word cName_;
            
            //- Is the coupling variable passive or not? (chemical source term)
            Switch passive_;

            //- when this variable should be mixed (two step mixing)?
            label mixStep_;

    public:
    
        //- Runtime type information
        TypeName("couplingVar");
            
        // Constructors
                       
            //- Construct from dictionary and mesh
            couplingVar
            (
                const dictionary& entry,
                const fvMesh& mesh
            );
            
            //- Construct copy
            couplingVar(const couplingVar& cm);


        //- Destructor
        virtual ~couplingVar()
        {};
            
        // Member functions

            //- Return variable name
            inline word cName() const
            {
                return cName_;
            }                 

            //- Return the const coupling variable field
            inline Switch passive() const
            {
                return passive_;
            }        

            //- Return the const coupling variable mixing step
            inline label mixStep() const
            {
                return mixStep_;
            }  

	    virtual autoPtr<basicMMCVar> clone() const
	    {
		return autoPtr<basicMMCVar>(new couplingVar(*this));
	    }            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
