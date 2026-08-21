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
    Foam::mmcVarSet
    
\*---------------------------------------------------------------------------*/

#ifndef mmcVarSet_H
#define mmcVarSet_H

#include "varSet.H"
#include "basicMMCVar.H"
#include "couplingVar.H"
#include "referenceVar.H"
#include "crVar.H"

namespace Foam
{

typedef HashTable<label, word> labelHashTable;

/*---------------------------------------------------------------------------*\
                            Class VarSet Declaration
\*---------------------------------------------------------------------------*/

class mmcVarSet
:
    public varSet    
{

    protected:            

        //- hash to return the index of reference variable in the full set Xi 
        //- from Reference Name           
        HashTable<label, word> rVarIndexinXi_;
        
        //- hash to return the index of coupling variable in the full set Xi 
        //- from Coupling Name          
        HashTable<label, word> cVarIndexinXi_;
        
        //- hash to return the index of reference variable in the in the subset 
        //- of reference variables using Reference Name           
        HashTable<label, word> rVarIndexinXiR_;
        
        //- hash to return the index of coupling variable in the in the subset 
        //- of coupling variables using Coupling Name          
        HashTable<label, word> cVarIndexinXiC_;
               
    public:

        TypeName("mmcVarSet");

        // Constructors
        
        
            //- Null constructor.
            inline mmcVarSet();
            
            //- Construct from dictionary and mesh
            mmcVarSet
            (
                const dictionary& mmcVarDict,
                const fvMesh&
            );
            
            //- Construct copy
            mmcVarSet(const mmcVarSet& mmcVs);


        //- Destructor
        virtual ~mmcVarSet() {};

        // Member functions

            //- Return a table with index of reference variables in the 
            inline const HashTable<label, word>& rVarInXi() const;
            
            //- Return the list of coupling variable names
            inline const HashTable<label, word>& cVarInXi() const;
            
            //- Return the list of reference variable names
            inline const HashTable<label, word>& rVarInXiR() const;
            
            //- Return the list of coupling variable names
            inline const HashTable<label, word>& cVarInXiC() const;
                        
            //- Return the number of reference variables
            inline label numXiR() const;
            
            //- Return the number of coupling variables
            inline label numXiC() const;
            
            //- Size function this allows to use forAll macros            
            inline label size() const;
            
        // Static Member Functions

            //- Return a null mmcVarSet
            inline static const mmcVarSet& null();

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mmcVarSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
