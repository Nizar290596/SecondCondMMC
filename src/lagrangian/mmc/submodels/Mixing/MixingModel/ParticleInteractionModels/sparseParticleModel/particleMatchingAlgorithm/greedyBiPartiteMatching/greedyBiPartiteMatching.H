/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::greedyBiPartiteMatchin
 
Description
    A greedy algorithm to generate bi-partite matching of particles by limiting
    the maximum mixing distance

SourceFiles
    greedyBiPartiteMatching.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/



#ifndef greedyBiPartiteMatching_H
#define greedyBiPartiteMatching_H

#include "dictionary.H"
#include "kdTree.H"
#include "undirectedGraph.H"
#include <unordered_map>
#include "Tuple2.H"
#include "particleMatchingAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class greedyBiPartiteMatching Declaration
\*---------------------------------------------------------------------------*/

template<class particleType>
class greedyBiPartiteMatching
: 
    public particleMatchingAlgorithm<particleType>
    //public particleMatchingAlgorithm<eulerianFieldData>
{
    private:

        //- Normalisation factors
        const scalar ri_;           // for physical coordinates
        const List<scalar> Xii_;    // for reference variable

        //- Maximum mixing distance
        scalar maxMixDist_;

        //- Unidirected graph
        undirectedGraph graph_;

        //- Construct a partially connected graph from a set of vertices
        void constructGraph(const List<List<scalar>>& vertices);


    public:

        //- Constructor 
        greedyBiPartiteMatching
        (
            const scalar& ri,
            const List<scalar>& Xii,
            const dictionary& dict
        );


        //- Find particle pairs from the given list of vertices
        //  The particleType class must have following public member functions
        //  vector& position() :  Returns the position in x,y,z coordinates
        //  List<scalar> XiR() :  Returns list of reference variables 
        void findPairs
        (
            const DynamicList<particleType>& eulerianFieldList,
            //const DynamicList<eulerianFieldData>& eulerianFieldList,
            DynamicList<List<label>>& pairs
        ) override;

};
}   // End namespace


#include "greedyBiPartiteMatching.C"

#endif
