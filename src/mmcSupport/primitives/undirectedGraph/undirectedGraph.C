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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "undirectedGraph.H"
#include "HashSet.H"

void Foam::undirectedGraph::addVertex
(
    const label v,              // vertex to add
    const List<label>& edges,    // vertices connected to this vertex v
    const List<scalar>& weights  // weight of each edge
)
{
    // Check if vertex is already present in the list
    auto it = graph_.find(v);
    if (it != graph_.end())
    {
        auto& knownEdges = it->second;
        // Check if the edges need to be updated
        for (label i=0; i < edges.size(); i++)
        {
            bool found = false;
            for (auto& e : knownEdges)
            {
                if (e.first == edges[i])
                    found = true;
            }
            if (!found && edges[i] != v)
            {
                knownEdges.emplace_back(edges[i],weights[i]);
                // Update connected vertex edge list as well
                List<label> tempEdges(1,v);
                List<scalar> tempWeights(1,weights[i]);
                this->addVertex(edges[i],tempEdges,tempWeights);
            }
        }
    }
    // If vertex is not yet in the graph
    else
    {
        // Create the list of edges with respetive weights
        std::vector<std::pair<label,scalar>> weightedEdges;
        for (label i=0; i < edges.size(); i++)
        {
            if (edges[i] != v)
                weightedEdges.emplace_back(edges[i],weights[i]);
        }

        graph_.emplace(v,weightedEdges);

        // Update the connected edge lists
        for (label i=0; i < edges.size(); i++)
        {
            List<label> tempEdges(1,v);
            List<scalar> tempWeights(1,weights[i]);
            this->addVertex(edges[i],tempEdges,tempWeights);
        }
    }
}


const std::vector<std::pair<Foam::label,Foam::scalar> >&
Foam::undirectedGraph::getNeighbor(const label v) const
{
    auto it = graph_.find(v);
    if (it == graph_.end())
    {
        FatalError << "Element "<<v<<" not found in graph" << exit(FatalError);
        return graph_.begin()->second;
    }

    return it->second;
}


Foam::DynamicList<Foam::undirectedGraph::edgeContainer>
Foam::undirectedGraph::getEdges()
{
    DynamicList<edgeContainer> edges;

    // Vertex pair is defined to be in ascending order
    HashSet<Pair<label>> foundEdges;

    for (auto& v : graph_)
    {
        auto& connectedVertices = v.second;
        for (auto& e :  connectedVertices)
        {
            Pair<label> vertexPair(v.first,e.first);
            if (v.first > e.first)
                vertexPair.flip();
            if (foundEdges.insert(vertexPair))
            {
                edgeContainer container;
                if (v.first > e.first)
                {
                    container.p = e.first;
                    container.q = v.first;
                }
                else
                {
                    container.p = v.first;
                    container.q = e.first;
                }

                container.w = e.second;
                edges.append(container);
            }
        }
    }

    // Sort the list by edge weight
    std::sort(edges.begin(),edges.end());


    return edges;
}



