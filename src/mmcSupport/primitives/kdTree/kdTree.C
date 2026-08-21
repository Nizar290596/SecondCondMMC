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

#include "kdTree.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * *

template<class VectorType, int bucketSize>
Foam::kdTree<VectorType,bucketSize>::kdTree
(
    const List<VectorType>& particles,
    List<scalar> weights,
    const bool medianBasedSorting,
    const bool normalize
)
:
    particles_
    (
        particles
    ),
    DIM_(particles_[0].size()),
    nParticles_(particles_.size()),
    medianBasedSorting_(medianBasedSorting),
    normalize_(normalize)
{
    // if default argument of weights is given use unity
    if (weights.size() == 0)
        weights.resize(DIM_,1.0);
    
    constructTree(weights);
}


template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::constructTree
(
    const List<scalar>& weights
)
{
    particleIndex_.resize(particles_.size());
    
    // fill particleIndex_ from 0 to n 
    std::iota(particleIndex_.begin(),particleIndex_.end(),0);

    // calculate the number of boxes and allocate 
    // storage for them
    label m=1;
    for (size_t nTmp = nParticles_; nTmp; nTmp /= 2)
        m *= 2; // multiply with 2

    label nBoxes = 2*nParticles_ - (m/2);
    if (m < nBoxes)
        nBoxes = m;

    boxes_.reserve(nBoxes);

    // create the coordinates array as a contineous list
    // of coordinates, first all x coordinates, than all y
    // up to the dimensions
    scalar* coords =  new scalar[DIM_*nParticles_];

    for (label dim=0, kk=0; dim < DIM_; dim++, kk+= nParticles_)
    {
        for (label k=0; k < nParticles_; k++)
            coords[kk+k] = this->particleList()[k][dim];
    }

    // find the minimum and maximum in each dimension to span up the root box
    Pair<ndPoint> loHi = getMinMaxExtend(weights);
    
    norm_ = weights;
    if (normalize_)
    {
        // set the normalisation factor
        for (label i=0; i < DIM_; ++i)
        {
            norm_[i] *= (loHi.second()[i]-loHi.first()[i]);
            // if the norm is very close to zero, this dimension probably does 
            // not exist in the data set. This means that for calcSplitDimension
            // this could be used for slight numerical inaccuracies due to the 
            // division of a value close to zero
            // --> Avoid selecting this dimension by setting a large norm value
            if (norm_[i] < SMALL)
                norm_[i] = GREAT;
        }
    }
    boxes_.emplace_back(loHi.first(),loHi.second(),-1,-1,-1,0,nParticles_);


    // instead of a recursive function call, stack the to
    // do tasks
    std::stack<label> boxToProcess;

    boxToProcess.push(0);

    while (!boxToProcess.empty())
    {
        // current parten box and dimension
        label parentBoxI = boxToProcess.top();
        
        boxToProcess.pop();
     
        // calculate dimension to split based on the spread in
        // each dimension
        label splitDim = calcSplitDimension(boxes_[parentBoxI]);
        boxes_[parentBoxI].splitDim = splitDim;
        
        // particle index low and high
        auto pIndexLo = boxes_[parentBoxI].particleIndexLo;
        auto pIndexHi = boxes_[parentBoxI].particleIndexHi;

        // sort particle indices according to dimension
        // returns index to particle mid point
        label pMidPointI = 0;
        if (medianBasedSorting_)
            pMidPointI = sortIndices
            (
                pIndexLo, 
                pIndexHi,
                &coords[splitDim*nParticles_]
            );
        else
        {
            // find average value of the split dimension for all 
            // particles in the box
            scalar sum = 0;
            for (int i=pIndexLo; i < pIndexHi; ++i)
                sum += coords[splitDim*nParticles_ + particleIndex_[i]];
            
            const scalar aveInSplitDim = sum/(pIndexHi-pIndexLo);
            
            pMidPointI = sortIndicesBasedOnThreshold
            (
                pIndexLo, 
                pIndexHi,
                &coords[splitDim*nParticles_],
                aveInSplitDim
            );
        }
        
        // create the two daughter boxes left and right
        ndPoint hi = boxes_[parentBoxI].hi;
        ndPoint lo = boxes_[parentBoxI].lo;
        hi[splitDim] = coords[splitDim*nParticles_+particleIndex_[pMidPointI]];
        lo[splitDim] = coords[splitDim*nParticles_+particleIndex_[pMidPointI]];


        // left box
        boxes_.emplace_back
        (
            boxes_[parentBoxI].lo,hi,
            parentBoxI,-1,-1,
            pIndexLo,pMidPointI
        );
       
        shrinkBoxToFit(boxes_[boxes_.size()-1]);

        // set the parent box daughter index
        boxes_[parentBoxI].left = boxes_.size()-1;

        if ((pMidPointI - pIndexLo) >= 2.0*bucketsize_)
            boxToProcess.push(boxes_[parentBoxI].left);
        

        // right box
        boxes_.emplace_back
        (
            lo, boxes_[parentBoxI].hi,
            parentBoxI,-1,-1,
            pMidPointI, pIndexHi
        );
        
        shrinkBoxToFit(boxes_[boxes_.size()-1]);

        // set the parent box daughter index
        boxes_[parentBoxI].right = boxes_.size()-1;
       
        // if size of elements between low and high is larger 
        // than bucket size split
        if ((pIndexHi - pMidPointI) >= 2.0*bucketsize_)
            boxToProcess.push(boxes_[parentBoxI].right);
    }

    delete[] coords;
}


template<class VectorType, int bucketSize>
Foam::label Foam::kdTree<VectorType,bucketSize>::calcSplitDimension(const Box& box)
{
    scalar maxSpread = 0;
    label spreadDim = -1;
    for (int i=0; i < DIM_; ++i)
    {
        const scalar spread = (box.hi[i]-box.lo[i])/norm_[i];
        if (spread > maxSpread)
        {
            maxSpread = spread;
            spreadDim = i;
        }
    }

    return spreadDim;
}


template<class VectorType, int bucketSize>
Foam::label
Foam::kdTree<VectorType,bucketSize>::sortIndicesBasedOnThreshold
(
    const label pIndexLo, 
    const label pIndexHi,
    const scalar* coords,
    const scalar threshold
)
{
    //-  Move indices in ind[li..ri] so that the elements in [li .. k] 
    //-  are less than the [k+1..ri] elements
   
    label k = pIndexLo; // left running index
    label ri = pIndexHi-1; // right running index

    while (k < ri)
    {
        if ((coords[particleIndex_[k]]) <= threshold)
            k++; //- good where it is.
        else
        {
            std::swap(particleIndex_[k],particleIndex_[ri]);
            ri--;
        }
    }

    if ((coords[particleIndex_[k]]) > threshold)
        k = k -1;

    // Mid point points to index of k+1
    return k+1;
}


template<class VectorType, int bucketSize>
Foam::label Foam::kdTree<VectorType,bucketSize>::sortIndices
(
    label pIndexLo,
    label pIndexHi,
    const scalar* coords
)
{
    // number of particles to consider
    const label n = pIndexHi - pIndexLo;
    
    // sort up to the mid point of the array
    label k = (pIndexHi - pIndexLo)/2;
    if ((pIndexHi - (k+pIndexLo)) % 2 != 0)
        k++;
    
    // start point of the particle indices 
    label* indx = &particleIndex_[pIndexLo];
    
    // temporary mid point variable
    label mid;

	label l=0;      // left element
	label ir=n-1;   // right element
	for (;;) 
    {
		if (ir <= l+1) 
        {
            // check between the last two elements
			if (ir == l+1 && coords[indx[ir]] < coords[indx[l]])
				std::swap(indx[l],indx[ir]);
            // return the index in particleIndex_
            // return the point after k+1
			return pIndexLo + k;
		} 
        else 
        {
			mid=(l+ir)/2;
			std::swap(indx[mid],indx[l+1]);
			if (coords[indx[l]] > coords[indx[ir]]) 
                std::swap(indx[l],indx[ir]);
			if (coords[indx[l+1]] > coords[indx[ir]]) 
                std::swap(indx[l+1],indx[ir]);
			if (coords[indx[l]] > coords[indx[l+1]]) 
                std::swap(indx[l],indx[l+1]);
			label i = l+1;
			label j = ir;
			label ia = indx[l+1];
			scalar a = coords[ia];
			for (;;) 
            {
				do i++; while (coords[indx[i]] < a);
				do j--; while (coords[indx[j]] > a);
				if (j < i) break;
				std::swap(indx[i],indx[j]);
			}
			indx[l+1]=indx[j];
			indx[j]=ia;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}




template<class VectorType, int bucketSize>
Foam::Pair<typename Foam::kdTree<VectorType,bucketSize>::ndPoint>
Foam::kdTree<VectorType,bucketSize>::getMinMaxExtend(const List<scalar>& weights)
{
    // Loop over particles to find min max value
    ndPoint minVal(DIM_,1E+15);
    ndPoint maxVal(DIM_,-1E+15);

    for (auto& p : particleList())
    {
        for (int i=0; i < DIM_; i++)
        {
            minVal[i] = std::min(minVal[i],p[i]);
            maxVal[i] = std::max(maxVal[i],p[i]);
        }
    }

    return Pair<ndPoint>(minVal,maxVal);
}


template<class VectorType, int bucketSize>
Foam::scalar Foam::kdTree<VectorType,bucketSize>::ndDist
(
    const VectorType& q,
    const VectorType& p
) const
{
    scalar dist = 0;

    for (int i=0; i < DIM_; ++i)
        dist += Foam::sqr((q[i]-p[i])/norm_[i]);

    return dist;
}


template<class VectorType, int bucketSize>
Foam::scalar 
Foam::kdTree<VectorType,bucketSize>::boxToPointDist
(
    const Box& box, 
    const VectorType& q
) const
{
    scalar d=0;
    for (int i =0; i < DIM_; i++)
    {
        if (q[i] < box.lo[i])
            d += Foam::sqr((q[i]-box.lo[i])/norm_[i]);

        if (q[i] > box.hi[i])
            d += Foam::sqr((q[i]-box.hi[i])/norm_[i]);
    }
     
    return d;
}


template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::csvWrite
(
    fileName name,
    label leftIndex,
    label rightIndex
)
{
    OFstream os(name);

    os << "# x,y,z,f"<<endl;
    // write query pos for kdTree
    for (int i = leftIndex; i < rightIndex; i++)
    {
        for (int k=0; k < DIM_-1; k++)
            os << particleList()[particleIndex_[i]][k] << ",";

        os << particleList()[particleIndex_[i]][DIM_-1] << endl;
    }
    os<<flush;
}


template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::shrinkBoxToFit(Box& box)
{
    ndPoint minVal(DIM_,1E+15);
    ndPoint maxVal(DIM_,-1E+15);
    for
    (
        label i=box.particleIndexLo;
        i < box.particleIndexHi;
        i++
    )
    {
        auto p = particleList()[particleIndex_[i]];
        for (int k=0; k < DIM_; ++k)
        {
            minVal[k] = std::min(minVal[k],p[k]);
            maxVal[k] = std::max(maxVal[k],p[k]);
        }
    }
    box.lo = minVal;
    box.hi = maxVal;
}


template<class VectorType, int bucketSize>
Foam::Tuple2<Foam::label,Foam::scalar> 
Foam::kdTree<VectorType,bucketSize>::nearestPtoP_inBox
(
    const Box& box,
    const label pI,
    const List<bool>& isPaired
) const 
{
    // ID of the found particle to pair
    label pairedParticle = -1;

    // Get the particle
    const VectorType& p = particleList()[pI];
    
    scalar minDist = GREAT;
    
    for
    (
        label i=box.particleIndexLo;
        i < box.particleIndexHi;
        i++
    )
    {
        if (particleIndex_[i] == pI || isPaired[particleIndex_[i]])
            continue;
        
        
        scalar dist = ndDist(p,particleList()[particleIndex_[i]]);
        if (dist < minDist)
        {
            minDist = dist;
            pairedParticle = particleIndex_[i];
        }
    }
    
    return Tuple2<label,scalar>(pairedParticle,minDist);
}


template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::searchInBox
(
    foundParticleQueue& queue,
    const VectorType& q,
    const Box& box,
    const label n
) const
{
    for
    (
        label i=box.particleIndexLo;
        i < box.particleIndexHi;
        i++
    )
    {
        scalar dist = ndDist(q,particleList()[particleIndex_[i]]);
        if (queue.size() < static_cast<size_t>(n))
        {
            queue.emplace
            (
                particleIndex_[i],
                dist
            );
        }
        else
        {
            if (dist < queue.top().dist)
            {
                queue.emplace
                (
                    particleIndex_[i],
                    dist
                );
            }
            if (queue.size() > static_cast<size_t>(n))
                queue.pop();
        }
    }
}


template<class VectorType, int bucketSize>
Foam::label Foam::kdTree<VectorType,bucketSize>::locate 
(
    const label pI
) const
{
    label boxIndexI =0;
    
    // Get the particle
    const VectorType& p = particleList()[pI];

    while 
    (
        boxes_[boxIndexI].left != -1 
     && boxes_[boxIndexI].right != -1 
    )
    {
        if (boxes_[boxIndexI].left == -1)
            boxIndexI = boxes_[boxIndexI].right;
        else if (boxes_[boxIndexI].right == -1)
            boxIndexI = boxes_[boxIndexI].left;
        else
        {
            label leftBoxI = boxes_[boxIndexI].left;
            label rightBoxI = boxes_[boxIndexI].right;

            label splitDim = boxes_[boxIndexI].splitDim;


            if (p[splitDim] <= boxes_[leftBoxI].hi[splitDim])
                boxIndexI = leftBoxI;
            else if (p[splitDim] > boxes_[rightBoxI].lo[splitDim])
                boxIndexI = rightBoxI;
            else
            {
                // check which one is closer
                scalar distLeft = p[splitDim] - boxes_[leftBoxI].hi[splitDim];
                scalar distRight = boxes_[rightBoxI].lo[splitDim] - p[splitDim];

                if (distLeft < distRight)
                    boxIndexI = leftBoxI;
                else
                    boxIndexI = rightBoxI;
            }
        }
    }
    
    return boxIndexI;
}

template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::findLeaf
(
    const Box& box, 
    std::deque<List<label>>& leafParticleList
) const
{
    // Check if it is a leaf node
    if (box.isLeaf())
    {
        DynamicList<label> pIndices;
        pIndices.reserve(box.particleIndexHi-box.particleIndexLo);
        for
        (
            label i=box.particleIndexLo;
            i < box.particleIndexHi;
            i++
        )
            pIndices.append(particleIndex_[i]);
        
        leafParticleList.push_back(pIndices);
    }
    else
    {
        if (box.left != -1)
            findLeaf(boxes_[box.left],leafParticleList);
        if (box.right != -1)
            findLeaf(boxes_[box.right],leafParticleList);
    }
}



// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * *

template<class VectorType, int bucketSize>
typename Foam::kdTree<VectorType,bucketSize>::resList
Foam::kdTree<VectorType,bucketSize>::nNearest
(
    const VectorType& q,
    const label n
) const
{
    // Loop over dimension and go through boxes
    
    //label splitDim = 0;
    label boxIndexI =0;
    

    while 
    (
        boxes_[boxIndexI].left != -1 
     && boxes_[boxIndexI].right != -1 
    )
    {
        if (boxes_[boxIndexI].left == -1)
            boxIndexI = boxes_[boxIndexI].right;
        else if (boxes_[boxIndexI].right == -1)
            boxIndexI = boxes_[boxIndexI].left;
        else
        {
            label leftBoxI = boxes_[boxIndexI].left;
            label rightBoxI = boxes_[boxIndexI].right;

            label splitDim = boxes_[boxIndexI].splitDim;


            if (q[splitDim] <= boxes_[leftBoxI].hi[splitDim])
                boxIndexI = leftBoxI;
            else if (q[splitDim] > boxes_[rightBoxI].lo[splitDim])
                boxIndexI = rightBoxI;
            else
            {
                // check which one is closer
                scalar distLeft = q[splitDim] - boxes_[leftBoxI].hi[splitDim];
                scalar distRight = boxes_[rightBoxI].lo[splitDim] - q[splitDim];

                if (distLeft < distRight)
                    boxIndexI = leftBoxI;
                else
                    boxIndexI = rightBoxI;
            }
        }
    }
    
    // reference to current parent box
    const Box& parentBox = boxes_[boxIndexI];

    // search in this box for nearest neighbors
    foundParticleQueue queue;

    if (parentBox.left != -1 || parentBox.right != -1)
        FatalError <<"Parent box is not right"<<exit(FatalError);

    searchInBox(queue,q,parentBox,n);

    // ============ Bottom Up Approach ================

    // Go from the current node up the tree
    List<bool> visitedBox(boxes_.size(),false);

    visitedBox[boxIndexI] = true;

    std::stack<label> boxToProcess;
    boxToProcess.push(boxes_[boxIndexI].parent);

    while (!boxToProcess.empty())
    {
        const label boxI = boxToProcess.top();
        boxToProcess.pop();
       
        const auto& curBox = boxes_[boxI];

        visitedBox[boxI] = true;

        // if the parent has not been visited yet push the parent
        if (curBox.parent != -1 && !visitedBox[curBox.parent])
            boxToProcess.push(curBox.parent);

        // if it is a leaf box process points
        if ( curBox.left == -1 && curBox.right == -1)
        { 
            if (queue.size() < static_cast<size_t>(n) || boxToPointDist(curBox,q) < queue.top().dist)
                searchInBox(queue,q,curBox,n);
        }
        else
        {
            if (queue.size() < static_cast<size_t>(n) || boxToPointDist(curBox,q) < queue.top().dist)
            {
                if (!visitedBox[curBox.left] && curBox.left != -1)
                    boxToProcess.push(curBox.left);
                if (!visitedBox[curBox.right] && curBox.right != -1)
                    boxToProcess.push(curBox.right);
            }
        }
    }

    // convert queue into vector and calculate distance in each dimension
    std::vector<foundParticle> foundParticleList;
    foundParticleList.reserve(queue.size());
    while (!queue.empty())
    {
        foundParticleList.push_back(std::move(queue.top()));
        queue.pop();
    }

    // calcualte distance
    for (auto& res : foundParticleList)
    {
        res.disReal.resize(DIM_);
        for (label i=0; i < DIM_; i++)
            res.disReal[i] = std::abs(particleList()[res.idx][i]-q[i]);
    }

    return foundParticleList;
}





template<class VectorType, int bucketSize>
void Foam::kdTree<VectorType,bucketSize>::printTree()
{
    label boxI = 0;
    // loop over all boxes
    for (Box& box : boxes_)
    {
        csvWrite
        (
            "box-"+Foam::name(boxI++)+".csv",
            box.particleIndexLo,
            box.particleIndexHi
        );
    }
}


template<class VectorType, int bucketSize>
Foam::DynamicList<typename Foam::kdTree<VectorType,bucketSize>::particlePair>
Foam::kdTree<VectorType,bucketSize>::findUniquePairs() const
{
    DynamicList<particlePair> particlePairs;
    particlePairs.reserve(floor(0.5*nParticles_));
    
    // A list to check if the particle is already paired
    List<bool> isPaired(nParticles_,false);
    
    
    // Shuffle list to avoid bias
    List<label> indices(particles_.size());
    std::iota(indices.begin(),indices.end(),0);
    std::random_shuffle(indices.begin(),indices.end());
    
    // Use always a number divisble by 2 for finding pairs to avoid expensive
    // lookup through the tree for the last particle with no partner
    
    const label numParticlePairs = floor(0.5*nParticles_);
    
    
    // Loop over all particles in the list and find their closes neighbor
    for (label i=0; i < numParticlePairs*2; ++i)
    {
        // index of the particle to look at
        const label pI = indices[i];
        
        if (isPaired[pI])
            continue;
        
        // Locate box of the particle
        const label boxIndexI = locate(pI);
        
        // find closest point in box (if exists)
        auto particleToPair = nearestPtoP_inBox(boxes_[boxIndexI],pI,isPaired);
        
        if (particleToPair.first() != -1)
        {
            isPaired[pI] = true;
            isPaired[particleToPair.first()] = true;
            particlePairs.append
            (
                particlePair(pI,particleToPair.first(),particleToPair.second())
            );
            continue;
        }
        
        // if it is -1 no particle to pair was found in the box
        // use bottom up approach to find the next fitting particle
        
        // Go from the current node up the tree
        List<bool> visitedBox(boxes_.size(),false);

        visitedBox[boxIndexI] = true;

        std::stack<label> boxToProcess;
        if (boxes_[boxIndexI].parent != -1)
            boxToProcess.push(boxes_[boxIndexI].parent);

        // has the pairing partner been found
        bool particleFound = false;
        
        while (!boxToProcess.empty() && !particleFound)
        {
            const label boxI = boxToProcess.top();
            boxToProcess.pop();
           
            const auto& curBox = boxes_[boxI];

            visitedBox[boxI] = true;

            // if the parent has not been visited yet push the parent
            if (curBox.parent != -1 && !visitedBox[curBox.parent])
                boxToProcess.push(curBox.parent);

            // if it is a leaf box process points
            if ( curBox.left == -1 && curBox.right == -1)
            { 
                auto particleToPair = nearestPtoP_inBox(curBox,pI,isPaired);
                if (particleToPair.first() != -1)
                {
                    isPaired[pI] = true;
                    isPaired[particleToPair.first()] = true;
                    particlePairs.append
                    (
                        particlePair(pI,particleToPair.first(),particleToPair.second())
                    );
                    particleFound = true;
                }
            }
            else
            {
                // Push both daughter boxes to the list to process
                if (!visitedBox[curBox.left] && curBox.left != -1)
                    boxToProcess.push(curBox.left);
                if (!visitedBox[curBox.right] && curBox.right != -1)
                    boxToProcess.push(curBox.right);
            }
        }
    }
    
    return particlePairs;
}


template<class VectorType, int bucketSize>
Foam::List<Foam::List<Foam::label>> 
Foam::kdTree<VectorType,bucketSize>::leafParticles() const
{
    std::deque<List<label>> leafParticleQueueList;
    findLeaf(boxes_[0],leafParticleQueueList);
    
    // Create an OpenFOAM list from std::deque
    List<List<label>> leafParticleList(leafParticleQueueList.size());
    
    label i = 0;
    for (List<label>& e : leafParticleQueueList)
    {
        leafParticleList[i++] = std::move(e);
    }
    return leafParticleList;
}

