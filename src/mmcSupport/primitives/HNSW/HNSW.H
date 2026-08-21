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
    Foam::HNSW
 
Description
    A fast hierarchical navigable small world (HNSW) algorithm to find 
    k closest neighbors to a given query point in n dimensions. This 
    class constructor is written to be used with the mmcFoam particle
    search in the KernelEstimation class, hence the use of dims and 
    weights. The graph and search functions work however for any dataType
    structure. 

    Using the constructor for mmc will create a euclidian distance
    function
    \f[
        \vec{d} = \sqrt{\sum_n{\left(q_i-p_i\right)^2}}
    \f]

    You can easily provide any distance function you want by
    overwriting the public member variable ditFunc_ with e.g a lambda
    expression:

    \verbatim
    distFunc_ = [](const dataType& q, const dataType& p) -> scalar
    {
        Write your own equation here...
    };
    \endverbatim

   Reference:
    \verbatim
        [1] Malkov, Y. A., & Yashunin, D. A. (2016). 
            Efficient and robust approximate nearest neighbor search using 
            Hierarchical Navigable Small World graphs. arXiv. 
            https://doi.org/10.48550/arXiv.1603.09320
    \endverbatim

SourceFiles
    HNSW.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#ifndef HNSW_H
#define HNSW_H

#include "List.H"
#include "OFstream.H"
#include <random>
#include <vector>
#include <unordered_set>
#include <queue>

namespace Foam
{

template<class dataType>
class hnsw
{
    private:

        // Overload std::pair to have descriptive functions to access
        class Pair
        :
            public std::pair<scalar,label>
        {
            public:
                scalar dist() const {return this->first;}
                label node() const {return this->second;}
           
                Pair() = default;

                Pair(std::pair<scalar,label> p)
                {
                    this->first = p.first;
                    this->second = p.second;
                }

                Pair(scalar dist, label i)
                {
                    first=dist;
                    second=i;
                }
        };

        struct compareByFirstDescending
        {
            constexpr bool operator()(Pair const &a, Pair const &b) const noexcept 
            {
                return a.first < b.first;
            }
        };
        
        struct compareByFirstAscending
        {
            constexpr bool operator()(Pair const &a, Pair const &b) const noexcept 
            {
                return a.first > b.first;
            }
        };

        // Descending means top element is furthest away
        typedef std::priority_queue<Pair,std::vector<Pair>,compareByFirstDescending> nodeConnectionsDescending;
        
        // Ascending means top element is closest 
        typedef std::priority_queue<Pair,std::vector<Pair>,compareByFirstAscending> nodeConnectionsAscending;

        // Fixed constants

            //- Max number of edges
            const label maxM_;

            //- Max number of edges in layer 0
            const label maxM0_;

            //- Normalization factor for level generation
            const scalar mult_;

            //- Number of connections during construction
            const label efConstruction_;

        // Private member variables 
            
            //- Store reference to data list
            List<dataType> data_;

            //- A list of visited nodes, required for searchLayer
            bool* visitedNodes_;
            
            //- Dimensions of the used elements to use for the distance function
            //  This is a special function required for the handling of
            //  the particles of mmcFoam
            const labelList dims_;

            //- List of weights to weigh the dimensions for the distance 
            //  function 
            const scalarList weights_;

            //- For each element in data_ a list of links to other elements
            //  Acess is: linkedLists[element+layer]
            //  First entry is distance, second entry is node label 
            List<std::vector<Pair>> linkedLists_;

            //- Random generator for level
            std::default_random_engine level_generator_;
            
            //- The entry point to the top layer HNSW graph
            label enter_point_=-1;

            //- Max level of the current graph
            label maxLevel_=-1;

        // Private Member Functions
            
            //- Add a new node to the graph with the label i of data_
            void addPoint(const label i, const label level);

            //- Search a layer for the ef closest points to a query point q
            //  Returns list of pairs of distances and labels 
            std::vector<Pair> searchLayer
            (
                const dataType& q, 
                const label enterPoint, 
                const label ef, 
                const label level
            );

            //- Select neighbors
            std::vector<Pair> selectNeighbors
            (
                const dataType& q,
                const std::vector<Pair>& C,
                const label M,
                const label level,
                bool extendCandidates=false,  // by default false
                bool keepPrunedConnections=true
            );

            //- Update the neighbor list of a point q
            void updateNeighbors
            (
                const label q,
                const label i,
                const label level
            );

            //- Remove the connection in the edges list to node b
            void removeConnection
            (
                std::vector<Pair>& edges,
                const label n
            );

            //- Return a random layer level
            label getRandomLevel();

            //- Write particle data as VTK legacy file format
            //  Useful for visualization and debugging
            void writeVTKLegacy
            (
                const fileName& name,
                const std::vector<label>& pointIDs,
                const std::vector<std::vector<label>>& edges
            );

    public:

        // Public member variable distFunc to allow access to it 
        std::function<scalar(const dataType&, const dataType&)> distFunc_;

        //- Constructor for OpenFOAM particle list 
        //  Requires special handling because data points and query points
        //  do not have the same structure
        hnsw
        (
            const List<dataType>& data, // For this dataType has to be a scalarList
            const labelList& dims,
            const scalarList& weights,
            const label M=16,          // A range between 5<M<48 is reasonable
            const label efConstruction=100
        );

        ~hnsw()
        {
            delete[] visitedNodes_;
        }

        //- Search n nearest neighbors to query point q
        labelList nNearestNeighbors(const dataType& q, const label n);

        //- Search n nearest neighbors to query point q and print path as vtk
        labelList nNearestNeighborsAndWriteVTK
        (
            const dataType& q,
            const label n,
            const fileName name = "nodePath.vtk"
        );

};
};  // End of namespace Foam

#include "HNSW.C"


#endif
