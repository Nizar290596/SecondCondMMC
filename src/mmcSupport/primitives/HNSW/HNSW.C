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

\*---------------------------------------------------------------------------*/

#include "HNSW.H"

template<class dataType>
void Foam::hnsw<dataType>::addPoint(const label i, const label level)
{
    // First check that label i is less than the size of data_
    if (i >= data_.size())
        FatalErrorInFunction() 
            << "Inserted point is beyond data size limit"
            << exit(FatalError);

    // Check if point exists by the number of connections in layer 0
    // if it is zero add the point
    if (linkedLists_[i].size() > 0)
        return;

    // Get layer to insert new data point
    label curLevel = getRandomLevel();

    if (level > 0)
        curLevel = level;

    label currNode = enter_point_;
    label currEnterPoint= enter_point_;

    if (currNode != -1)
    {
        if (curLevel < maxLevel_)
        {
            for (label level = maxLevel_; level > curLevel; level--)
            {    
                std::vector<Pair> W = searchLayer
                (
                    data_[i],
                    currEnterPoint,
                    1,
                    level
                );
                // W is in descending order meaning the first element if the 
                // furthest away
                currEnterPoint = W.back().node();                
            }
        }

        for (label level = min(curLevel,maxLevel_); level >= 0; --level)
        {
            auto W = searchLayer(data_[i],currEnterPoint,efConstruction_,level);
            const label M = (level == 0) ? maxM0_ : maxM_;
            auto neighbors = selectNeighbors(data_[i],W,M,level,false,true);

            // Add connections between neighbors and q aka data_[i]
            for (auto e : neighbors)
            {
                // If adding the link would leed to more edges than
                // allowed --> Shrink list
                if (linkedLists_[e.node()+level].size() >= M)
                {
                    updateNeighbors(e.node(),i,level);
                }
                else
                {
                    linkedLists_[e.node()+level].emplace_back
                    (
                        e.dist(),
                        i
                    );

                    linkedLists_[i+level].push_back(e);
                }
            }
        }
    }
    else
    {
        enter_point_ = i;
        maxLevel_ = curLevel;
    }

    if (curLevel > maxLevel_)
    {
        enter_point_ = i;
        maxLevel_ = curLevel;
    }
}
        

template<class dataType>
std::vector<typename Foam::hnsw<dataType>::Pair>
Foam::hnsw<dataType>::searchLayer
(
    const dataType& q,
    const label enterPoint,
    const label ef,
    const label level
)
{
    // Set enter point of visited nodes to true
    visitedNodes_[enterPoint] = true;

    // List of candidates in a ascending order
    // top() returns the element closest to q;
    nodeConnectionsAscending C;
    C.emplace
    (
        distFunc_(data_[enterPoint],q),
        enterPoint
    );

    // List of found nearest neighbours in a descending order
    // top() returns the element furthest to q;
    nodeConnectionsDescending W;
    W.emplace
    (
        distFunc_(data_[enterPoint],q),
        enterPoint
    );

    while (C.size() > 0)
    {
        // extract nearest element to q in C
        auto c = C.top();
        C.pop();
        scalar maxDistCq = c.dist();

        // get furthest element to q in W
        scalar maxDistWq = W.top().dist();

        // if the distance in candidates is larger
        // than in W --> break while loop
        if (maxDistCq > maxDistWq)
            break;
        
        // Loop over all edges of C
        for (auto& e : linkedLists_[c.node()+level])
        {
            if (!visitedNodes_[e.node()])
            {
                visitedNodes_[e.node()]=true;

                // get furthest element to q in W
                scalar maxDistWq = W.top().dist();
                scalar dist = distFunc_(data_[e.node()],q);
                
                if (dist < maxDistWq)
                {
                    C.emplace
                    (
                        dist,
                        e.node()
                    );

                    W.emplace
                    (
                        dist,
                        e.node()
                    );
                    
                    if ( W.size() > ef )
                    {
                        W.pop();
                    }
                }
            }
        }
    }

    std::vector<Pair> Wvec;
    Wvec.reserve(W.size());
    while (!W.empty())
    {
        Wvec.emplace_back(std::move(W.top()));
        W.pop();
    }

    // Clear visited nodes
    memset(visitedNodes_, 0, sizeof(bool) * data_.size());

    return Wvec;
}


template<class dataType>
typename std::vector<typename Foam::hnsw<dataType>::Pair>
Foam::hnsw<dataType>::selectNeighbors
(
    const dataType& q,
    const std::vector<Pair>& C,
    const label M,
    const label level,
    bool extendCandidates,  // by default false
    bool keepPrunedConnections
)
{
    // if extendCandidates is turned off and 
    // keepPrunedConnections is turned on the vector R
    // if effectively filled by the M closest values of C
    std::vector<Pair> R;
    R.reserve(M);

    // working queue for the candidates 
    // top() returns closest 
    nodeConnectionsAscending W(C.begin(),C.end());

    if (!extendCandidates && keepPrunedConnections)
    {
        while (W.size() > 0 && R.size() < M)
        {
            auto c = W.top();
            W.pop();
            R.push_back(std::move(c));
        }
    }
    else
    {
        // Expand candidate list
        if (extendCandidates)
        {
            for (auto& e : C)
            {
                for (auto& eAdj : linkedLists_[e.node()+level])
                {
                    W.emplace
                    (
                        distFunc_(q,data_[e.node()]),
                        eAdj.node()
                    );
                }
            }
        }

        
        // queue of discarded candidates
        // top() returns closest 
        nodeConnectionsAscending Wd;

        // Minimum distance in R to q
        scalar minDistRQ = GREAT;

        while (W.size() > 0 && R.size() < M)
        {
            auto c = W.top();
            W.pop();

            if (c.dist() < minDistRQ)
            {
                R.push_back(c);
                minDistRQ = c.dist();
            }
            else
            {
                Wd.push(c);
            }
        }

        if (keepPrunedConnections)
        {
            while (Wd.size() > 0 && R.size() < M)
            {
                auto c = Wd.top();
                Wd.pop();
                R.push_back(std::move(c));
            }
        }
    }

    return R;
}


template<class dataType>
void Foam::hnsw<dataType>::updateNeighbors
(
    const label q,
    const label i,
    const label level
)
{
    auto& edges = linkedLists_[q+level];

    // Distance from node q to i
    const scalar dist0 = distFunc_(data_[q],data_[i]);
 
    std::sort(edges.begin(),edges.end());

    if (edges.back().dist() > dist0)
    {
        const label oldNeighbor = edges.back().node();

        // Remove connection in oldNeighbor
        removeConnection
        (
            linkedLists_[oldNeighbor+level],
            q
        );

        // Add connection to list for point q
        edges.emplace_back
        (
            dist0,
            i
        );

        // Add connection to list for point i
        linkedLists_[i+level].emplace_back
        (
            dist0,
            q
        );
    }
}


template<class dataType>
Foam::label Foam::hnsw<dataType>::getRandomLevel()
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double r = -log(distribution(level_generator_)) * mult_;
    return label(r);
}


template<class dataType>
void Foam::hnsw<dataType>::removeConnection
(
    std::vector<Pair>& edges,
    const label n
)
{
    for (auto it = edges.begin(); it != edges.end(); )
    {
        // Remove entry if node exists
        if (it->node() == n)
            it = edges.erase(it);
        else
            ++it;
    }
}


template<class dataType>
void Foam::hnsw<dataType>::writeVTKLegacy
(
    const fileName& name,
    const std::vector<label>& pointIDs,
    const std::vector<std::vector<label>>& edges
)
{

    OFstream os(name);
    
    // Write header
    os << "# vtk DataFile Version 1.0" << nl
       << "node connections"<< nl
       << "ASCII" << nl << nl
       << "DATASET POLYDATA" << nl
       << endl;

    // Write points
    os << "POINTS "<<pointIDs.size()<<" float"<<endl;
    for (auto id : pointIDs)
    {
        os << data_[id][0] << " "
           << data_[id][1] << " "
           << data_[id][2] << endl;
    }
    os << nl;

    // calculate all values of the lines list
    label nLineEntries = edges.size();
    forAll(edges,i)
    {
        nLineEntries += edges[i].size();
    }

    // Write connections
    os << "LINES " << edges.size() << " " << nLineEntries<<endl;
    for (auto& e : edges)
    {
        os << e.size() << " ";
        for (int i=0; i < e.size()-1; ++i)
        {
            os << e[i] << " ";
        }
        os << e[e.size()-1] << endl;
    }
}
// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * 

template<class dataType>
Foam::hnsw<dataType>::hnsw
(
    const List<dataType>& data, // For this dataType has to be a scalarList
    const labelList& dims,
    const scalarList& weights,
    const label M,               // A range between 5<M<48 is reasonable
    const label efConstruction
)
:
maxM_(M),   // From reference 
maxM0_(2*M), // From reference
mult_(1.0/std::log(1.0*M)),
efConstruction_(std::max(efConstruction,M)),
dims_(dims),
weights_(weights)
{
    // Use 100 as a random seed for the level generator
    level_generator_.seed(100);

    visitedNodes_ = new bool[data.size()];

    linkedLists_.resize(data.size()*label((1.0/mult_)));
    data_.resize(data.size());

    forAll(data,i)
    {
        scalarList element(dims.size());
        forAll(dims,j)
        {
            element[j] = weights[j]*data[i][dims[j]];
        }
        data_[i] = std::move(element);
    }

    distFunc_ = [](const dataType& e1, const dataType& e2) -> scalar
    {
        scalar dist = 0;
        for (int i=0; i<e1.size(); ++i)
        {
            dist += (e1[i]-e2[i])*(e1[i]-e2[i]);
        }
        return dist;
    };

    // Add all data points to graph
    forAll(data_,i)
    {
        addPoint(i,-1);
    }
}


template<class dataType>
Foam::labelList Foam::hnsw<dataType>::nNearestNeighbors
(
    const dataType& q,
    const label n
)
{
    label currEnterPoint = enter_point_;
    for (label level=maxLevel_; level > 0; --level)
    {
        auto W = searchLayer(q,currEnterPoint,1,level);
        // W is in descending order and first element
        // is furthest away
        currEnterPoint=W.back().node();
    }

    auto W = searchLayer(q,currEnterPoint,n,0);

    labelList neighbors(W.size());
    
    forAll(neighbors,i)
    {
        neighbors[i] = W[i].node();
    }


    return neighbors;
}


template<class dataType>
Foam::labelList Foam::hnsw<dataType>::nNearestNeighborsAndWriteVTK
(
    const dataType& q,
    const label n,
    const fileName name
)
{
    std::vector<label> pointIDs;
    std::vector<std::vector<label>> edges;

    label currEnterPoint = enter_point_;

    pointIDs.push_back(currEnterPoint);

    for (label level=maxLevel_; level > 0; --level)
    {
        auto W = searchLayer(q,currEnterPoint,1,level);
        // W is in descending order and first element
        // is furthest away
        currEnterPoint=W.back().node();
        pointIDs.push_back(currEnterPoint);
    }

    // This search path is the first edge
    std::vector<label> edge;
    for (int i=0; i < pointIDs.size(); ++i)
        edge.push_back(i);

    edges.push_back(std::move(edge));

    // Save the position in pointIDs of the current enter point
    label enterPointInd = pointIDs.size()-1;

    auto W = searchLayer(q,currEnterPoint,n,0);

    labelList neighbors(W.size());
    
    forAll(neighbors,i)
    {
        neighbors[i] = W[i].node();
        pointIDs.push_back(W[i].node());
        std::vector<label> edge(2);
        edge[0] = enterPointInd;
        edge[1] = pointIDs.size()-1;
        edges.push_back(std::move(edge));
    }

    writeVTKLegacy
    (
        name,
        pointIDs,
        edges
    );

    return neighbors;
}

