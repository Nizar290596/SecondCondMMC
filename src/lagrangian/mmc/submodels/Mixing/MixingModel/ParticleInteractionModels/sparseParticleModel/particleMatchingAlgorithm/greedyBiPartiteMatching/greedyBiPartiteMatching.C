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
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY s
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "HashSet.H"

template<class particleType>
Foam::greedyBiPartiteMatching<particleType>::greedyBiPartiteMatching
(
    const scalar& ri,
    const List<scalar>& Xii,
    const dictionary& dict
)
: 
    particleMatchingAlgorithm<particleType>(ri,Xii,dict),
    ri_(ri),
    Xii_(Xii),
    maxMixDist_(dict.get<scalar>("maxMixingLength"))
{}

template<class particleType>
void Foam::greedyBiPartiteMatching<particleType>::constructGraph
(
    const List<List<scalar>>& vertices
)
{
    // Generate the weights
    List<scalar> weights(vertices.first().size(),1.0);
 
    // Generate the k-d tree to find nearest neighbors
    kdTree<List<scalar>> tree(vertices,weights,false,false);
    // Find for each vertex the m closest neighbors
    for (label i=0; i < vertices.size(); i++)
    {
        auto resList = tree.nNearest(vertices[i],80);
        // Convert res list to a format readable by the undirectedGraph
        DynamicList<label,50> connectedVertices;
        DynamicList<scalar,50> edgeWeights;
        for (auto& res : resList)
        {
            if (res.idx == i)
                continue;
            connectedVertices.append(res.idx);
            edgeWeights.append(res.dist);
        }

        graph_.addVertex(i,connectedVertices,edgeWeights);
    }
}



template<class particleType>
void Foam::greedyBiPartiteMatching<particleType>::findPairs
//void Foam::greedyBiPartiteMatching<eulerianFieldData>::findPairs
(
    const DynamicList<particleType>& eulerianFieldList,
    //const DynamicList<eulerianFieldData>& eulerianFieldList,
    DynamicList<List<label>>& pairs
)
{
    // Create a list of vertices from the list of eulerianDataFields
    // provided by mixParticleModel

    List<List<scalar>> vertices(eulerianFieldList.size());

    forAll(eulerianFieldList,i)
    {
        auto& e = eulerianFieldList[i];
        vertices[i].resize(3+e.XiR().size());
        vertices[i][0] = e.position()[0]/ri_;
        vertices[i][1] = e.position()[1]/ri_;
        vertices[i][2] = e.position()[2]/ri_;

        label k = 0;
        for (auto& refVar : e.XiR())
        {
            vertices[i][3 + k] = refVar/Xii_[k];
            k++;
        }
    }

    constructGraph(vertices);

    // get edge list
    auto edgeList = graph_.getEdges();

    // Keep track of particles that are already paired to avoid double
    // pairing
    HashSet<label> foundParticles;
    pairs.clear();
    for (auto& e: edgeList)
    {
        // if not paired yet, add to matching
        if (!foundParticles.test(e.p) && !foundParticles.test(e.q))
        {
            List<label> pair(2);
            pair[0] = e.p;
            pair[1] = e.q;
            pairs.append(pair);
            foundParticles.insert(e.p);
            foundParticles.insert(e.q);
        }
    }
}


