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
    Foam::referenceVar

Description
    This class is aimed to help the creation of reference in MMC. It may 
    also be used for variables required in postprocessing.

SourceFiles
    referenceVar.C

\*---------------------------------------------------------------------------*/

#ifndef referenceVar_H
#define referenceVar_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"

#include "basicMMCVar.H"
#include "evolModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class referenceVar Declaration
\*---------------------------------------------------------------------------*/

class referenceVar
:
    public virtual basicMMCVar
{
    protected:
              
        // Protected data
        
            //- coupling name (name in particles)
            word rName_;
            
            //- Is the coupling variable passive or not? (chemical source term)
            word refType_;

	    //based on refType build the evolution method
	    autoPtr<evolModel> evolMethod_; 
    public:
    
        //- Runtime type information
        TypeName("referenceVar");
        
        // Constructors
                   
            //- Construct from dictionary and mesh
            referenceVar
            (
                const dictionary& entry,
                const fvMesh& mesh
            );
            
            //- Construct copy
            referenceVar(const referenceVar& cm);

	    referenceVar& operator=(const referenceVar&) = delete; //unsure


        //- Destructor
        virtual ~referenceVar()
        {};
        
        // Member functions

            //- Return variable name
            inline word rName() const          
            {
                return rName_;
            }
            
            //- Return the const coupling variable field
            inline word refType() const
            {
                return refType_;
            }                          
            
	    //return reference variable evolution model
	    inline const evolModel& evolveMethod() const
	    {
		return evolMethod_();
	    }

	    //construct and return a clone to basic mmc variable
	    virtual autoPtr<basicMMCVar> clone() const
	    {
		return autoPtr<basicMMCVar>(new referenceVar(*this));
	    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
