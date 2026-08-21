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
    Foam::SinglePhaseMixture

Description
    Templated parcel single phase, multi-component class

SourceFiles
    SinglePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef SinglePhaseMixture_H
#define SinglePhaseMixture_H

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SinglePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SinglePhaseMixture
:
    public CompositionModel<CloudType>
{
    // Private data
    const Switch printMolFracEnabled_;
    // Private Member Functions

public:

    //- Runtime type information
    TypeName("singlePhaseMixture");


    // Constructors

        //- Construct from dictionary
        SinglePhaseMixture(const dictionary& dict, CloudType& owner);

    //- Destructor
    virtual ~SinglePhaseMixture();


    // Member Functions

    // Access

    // Mixture properties

        //- Return the list of component names for phaseI
        virtual const wordList& componentNames() const;

        //- Return the list of mixture mass fractions
        virtual const scalarField YMixture0(label patch,label patchFace);
        virtual const scalarField YMixture0(label celli);
        
        //- Return if mole fractions should be writen 
        virtual const Switch printMoleFractionsEnabled() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SinglePhaseMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
