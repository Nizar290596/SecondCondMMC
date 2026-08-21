/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::particleChemistryReductionMethods::pDRG

Description
    Implementation of the Directed Relation Graph (oDRG) method 
    for particle chemistry

SourceFiles
    pDRG.C

\*---------------------------------------------------------------------------*/

#ifndef pDRG_H
#define pDRG_H

#include "particleChemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace particleChemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                            Class pDRG Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class pDRG
:
    public particleChemistryReductionMethod<ReactionThermo, ThermoType>
{
    // Private data

        //- List of label for the search initiating set
        labelList searchInitSet_;

public:

    //- Runtime type information
    TypeName("pDRG");


    // Constructors

        //- Construct from components
        pDRG
        (
            const IOdictionary& dict,
            particleTDACChemistryModel<ReactionThermo, ThermoType>& chemistry
        );


    // Destructor
    virtual ~pDRG();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleChemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pDRG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
