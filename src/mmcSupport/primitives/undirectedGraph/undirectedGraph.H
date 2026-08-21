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
    Represent an undirected graph as a adjacency list. Store weight or distance
    of each edge

SourceFiles
    undirectedGraph.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/



#ifndef undirectedGraph_H
#define undirectedGraph_H

#include "dictionary.H"
#include <unordered_map>
#include "Tuple2.H"
#include "Pair.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class undirectedGraph Declaration
\*---------------------------------------------------------------------------*/

class undirectedGraph
{
    private:

        //- Store the graph as an adjacency list with weighted edges
        std::unordered_map
        <
            label,
            std::vector<std::pair<label,scalar> > 
        > graph_;

    public:
        
        //- Structure to store edge of vertex (p,q) and weight w 
        struct edgeContainer
        {
            label  p;
            label  q;
            scalar w;

            bool operator < (const edgeContainer& c2) const
            {
                return (w < c2.w);
            }
        };


    public:

        //- Empty default constructor
        undirectedGraph() = default;


        //- Add a vertex to the graph with a list of connections
        //  vertices in the edge list that are not yet part of the graph
        //  are added to the graph
        void addVertex
        (
            const label v,              // vertex to add
            const List<label>& edges,   // vertices connected to this vertex v
            const List<scalar>& weights // weight of each edge
        );

        //- Get all neighboring vertices to vertex v with their respective
        //  edge weights
        const std::vector<std::pair<label,scalar> >&
        getNeighbor(const label v) const;

        //- Return a list of edges in ascending order
        DynamicList<edgeContainer> getEdges();


    // Inline Functions

        inline label size() const {return graph_.size();}




};


}


#endif
