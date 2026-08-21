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
    Foam::basicMMCVar

Description
    This class is the base/interface handle reference and coupling 
    variables in MMC. It may olso be used for variables required in 
    postprocessing.

SourceFiles
    basicMMCVar.C

\*---------------------------------------------------------------------------*/

#ifndef basicMMCVar_H
#define basicMMCVar_H

#include "volFields.H"
#include "PtrList.H"
#include "hashedWordList.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "evolModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicMMCVar Declaration
\*---------------------------------------------------------------------------*/

class basicMMCVar
{
    protected:
           
        // Protected data
        
            //- LES name
            const word name_;

            //- coupling variable volumetric field
            autoPtr<volScalarField> field_;

    public:
    
        //- Runtime type information
        TypeName("basicMMCVar");
        
        //- Declare runtime constructor selection table    
        declareRunTimeSelectionTable
        (
            autoPtr,
            basicMMCVar,
            dictionary,
            (
                const dictionary& entry,
                const fvMesh& mesh
            ),
            (entry, mesh)
        );

        // Constructors
                   
            //- Construct from dictionary and mesh
            basicMMCVar
            (
                const dictionary& entry,
                const fvMesh& mesh
            );

	basicMMCVar(const basicMMCVar& bmv);

        //- Destructor
        virtual ~basicMMCVar()
        {};

        //- Selector
        static autoPtr<basicMMCVar> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );
        
        // Member functions

            //- Return variable name
            inline const word& name() const;

            //- Return the coupling variable field
            inline volScalarField& field();
            
            //- Return the const coupling variable field
            inline const volScalarField& field() const;
            
            //- virtual functions
            
            virtual inline word cName() const;
            
            virtual inline Switch passive() const;
            
            virtual inline label mixStep() const;
            
            virtual inline word rName() const;
            
            virtual inline word refType() const;

	    virtual inline const evolModel& evolveMethod() const;
            
            //- Construct and return a clone
            virtual autoPtr<basicMMCVar> clone() const
            {
                return autoPtr<basicMMCVar>(new basicMMCVar(*this));
            }

            
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicMMCVarI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
