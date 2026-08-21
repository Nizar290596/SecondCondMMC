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
    Foam::particleMatchingAlgorithm

Description
    Templated base class to select the particle matching algorithm
    Based on the idea of sub-modules 

SourceFiles
    particleMatchingAlgorithm.C
    particleMatchingAlgorithmNew.C

\*---------------------------------------------------------------------------*/

#ifndef particleMatchingAlgorithm_H
#define particleMatchingAlgorithm_H

#include "IOdictionary.H"
#include "autoPtr.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class particleMatchingAlgorithm Declaration
\*---------------------------------------------------------------------------*/

template<class particleType>
class particleMatchingAlgorithm
{
    public:

        // Constructors

            //- Construct from dictionary
            particleMatchingAlgorithm
            (
                const scalar& ri, 
                const List<scalar>& Xii,
                const dictionary& dict
            );


        //- Destructor
        virtual ~particleMatchingAlgorithm();


        //- Selector
        static autoPtr<particleMatchingAlgorithm> New
        (
            const scalar& ri, 
            const List<scalar>& Xii,
            const dictionary& dict
        );

        //- Find particle pairs from the given list of vertices
        //  The particleType class must have following public member functions
        //  vector& position() :  Returns the position in x,y,z coordinates
        //  List<scalar> XiR() :  Returns list of reference variables 
        virtual void findPairs
        (
            const DynamicList<particleType>& eulerianFieldList,
            //const DynamicList<eulerianFieldData>& eulerianFieldList,
            DynamicList<List<label>>& pairs
        ) 
        {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "particleMatchingAlgorithm.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
