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
    Foam::particleChemistryReductionMethods::pEFA

Description

SourceFiles
    pEFA.C

\*---------------------------------------------------------------------------*/

#ifndef pEFA_H
#define pEFA_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace particleChemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class pEFA
:
    public particleChemistryReductionMethod<ReactionThermo, ThermoType>
{
    // Private data

        List<label> sC_,sH_,sO_,sN_;
        scalar  sortPart_;

public:

    //- Runtime type information
    TypeName("pEFA");


    // Constructors

        //- Construct from components
        pEFA
        (
            const IOdictionary& dict,
            particleTDACChemistryModel<ReactionThermo, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~pEFA();


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
    #include "pEFA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
