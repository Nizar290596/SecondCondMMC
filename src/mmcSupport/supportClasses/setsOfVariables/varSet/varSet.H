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
    Foam::varSet
Description
    
SourceFiles
    varSet.C
    varSetI.H
    varSetIO.C
\*---------------------------------------------------------------------------*/

#ifndef varSet_H
#define varSet_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "basicMMCVar.H"
#include "couplingVar.H"
#include "referenceVar.H"
#include "crVar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef hashedWordList varTable;
/*---------------------------------------------------------------------------*\
                            Class varSet Declaration
\*---------------------------------------------------------------------------*/

class varSet
{

    protected:

        //- Table of variable names
        varTable varNames_;

        //- variable volumetric fields
        PtrList<basicMMCVar> Vars_;
        
    public:

        TypeName("varSet");

        // Constructors
        
            //- Null constructor.
            varSet();
            
            //- Construct from dictionary and mesh
            varSet
            (
                const dictionary& VarDict,
                const fvMesh&
            );
            
            //- Construct copy
            varSet(const varSet& vs);


        //- Destructor
        virtual ~varSet()
        {};

        // Member functions

            //- Return the list of variable names
            inline const varTable& varNames() const;

            //- Return the variable prt list
            inline PtrList<basicMMCVar>& Vars();

            //- Return the const variable ptr list
            inline const PtrList<basicMMCVar>& Vars() const;

            //- Return the a specific variable given by index
            inline basicMMCVar& Vars(const label i);

            //- Return the const field for a specific variable given by index
            inline const basicMMCVar& Vars(const label i) const;

            //- Return the field for a specific variable given by name
            inline basicMMCVar& Vars(const word& varName);

            //- Return the const field for a specific variable given by name
            inline const basicMMCVar& Vars(const word& varName) const;

            //- Does the list include this variable?
            inline bool contains(const word& varName) const;
            
        // Static Member Functions

            //- Return a null varSet
            inline static const varSet& null();

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "varSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
