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
    Foam::flameSheetParticleReaction

Description

SourceFiles
    NoReaction.C

\*---------------------------------------------------------------------------*/

#ifndef flameSheetParticleReaction_H
#define flameSheetParticleReaction_H

#include "ReactionModel.H"

#include "List.H"
#include "scalar.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class flameSheetParticleReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class flameSheetParticleReaction
:
    public ReactionModel<CloudType>
{
    // Private data

        //- Mixture fraction - stoichiometry, lower flammability limit, upper flammability limit
        const List<scalar> fSLU_;
        
        //- Quench temperature
        const scalar Tquench_;
        
        //-  Species mass fractions - stoichiometry, at f=0 and f=1
        List<scalarField> Y0S1_;

    public:

        //- Runtime type information
        TypeName("flameSheetParticleChemistry");


    // Constructors

        //- Construct from dictionary
        flameSheetParticleReaction(const dictionary& dict, CloudType& owner);

    //- Destructor
    virtual ~flameSheetParticleReaction();

    // Member functions
    
        //- Ignite the flamesheet
        virtual void ignite();
        
        //- Calcuate chemistry
        virtual void calculate
        (
            const scalar t0, 
            const scalar dt,
            const scalar ha, 
            const scalar pc, 
            const scalar pT,
            const scalar pz,
            scalarField& Y
        ) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "flameSheetParticleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
