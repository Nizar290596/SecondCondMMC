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

#include <algorithm>

template<class particleType>
Foam::kdTreeLikeMatching<particleType>::kdTreeLikeMatching
(
    const scalar& ri,
    const List<scalar>& Xii,
    const dictionary& dict
)
:
    particleMatchingAlgorithm<particleType>(ri,Xii,dict),
    //particleMatchingAlgorithm<eulerianFieldData>(ri,Xii,dict),
    ri_(ri),
    Xii_(Xii)
{}


template<class particleType>
void Foam::kdTreeLikeMatching<particleType>::findPairs
(
    const DynamicList<particleType>& eulerianFieldList,
    //const DynamicList<eulerianFieldData>& eulerianFieldList,
    DynamicList<List<label>>& pairs
)
{
    // Clear particle pairs first
    pairs.clear();

    // Keeping track of indices for kdTreeLikeSearch
    std::vector<label> L;
    std::vector<label> U;
    L.reserve(eulerianFieldList.size());
    U.reserve(eulerianFieldList.size());

    // create an index list for the particle data
    std::vector<label> pInd(eulerianFieldList.size());
    std::iota(pInd.begin(),pInd.end(),0);

    kdTreeLikeSearch(eulerianFieldList,1,eulerianFieldList.size(),pInd,L,U);

    // reserve space for list of pairs
    pairs.reserve(ceil(0.5*eulerianFieldList.size()));

    for(size_t i=0; i<L.size(); i++)
    {
        label p = L[i] - 1;

        label q = L[i];

        if(U[i] - L[i] < 2)
        {
            List<label> pair(2);
            pair[0] = pInd[p];
            pair[1] = pInd[q];

            pairs.append(std::move(pair));
        }
        else if(U[i] - L[i] == 2)
        {
            label r = L[i] + 1;

            List<label> pair(3);
            pair[0] = pInd[p];
            pair[1] = pInd[q];
            pair[2] = pInd[r];

            pairs.append(std::move(pair));
        }
    }
}


template<class particleType>
void Foam::kdTreeLikeMatching<particleType>::kdTreeLikeSearch
(
    const DynamicList<particleType>& particleList,
    //const DynamicList<eulerianFieldData>& particleList,
    label l,
    label u,
    std::vector<label>& pInd,
    std::vector<label>& L,
    std::vector<label>& U
) const
{
    //- Break the division if the particle list has length less than 2
    if (u - l <= 2)
    {
        //- Divide particles into groups of two or three
        L.push_back(l);

        U.push_back(u);

        return ;
    }

    label m = (l + u)/2;
    if ( (u - m) % 2 != 0 ) m++;

    auto iterL = pInd.begin();

    auto iterM = pInd.begin();

    auto iterU = pInd.begin();

    std::advance(iterL,l-1);

    std::advance(iterM,m-1);

    std::advance(iterU,u  );


    scalar maxInX = -GREAT;
    scalar maxInY = -GREAT;
    scalar maxInZ = -GREAT;
    List<scalar> maxInXiR(Xii_.size(),0.0);

    scalar minInX = GREAT;
    scalar minInY = GREAT;
    scalar minInZ = GREAT;
    List<scalar> minInXiR(Xii_.size(),GREAT);

    // Find minimum and maximum for each coordinate
    for (auto it = iterL; it != iterU; it++)
    {
        auto& pos = particleList[*it].position();
        maxInX = std::max(maxInX,pos.x());
        maxInY = std::max(maxInY,pos.y());
        maxInZ = std::max(maxInZ,pos.z());

        minInX = std::min(minInX,pos.x());
        minInY = std::min(minInY,pos.y());
        minInZ = std::min(minInZ,pos.z());

        forAll(Xii_,i)
        {
            maxInXiR[i] = std::max(maxInXiR[i],particleList[*it].XiR()[i]);
            minInXiR[i] = std::min(minInXiR[i],particleList[*it].XiR()[i]);
        }
    }

    //- Scaled/stretched distances between Max and Min in each direction
    //- Default is random mixing, overwritten if mixing distances greater than ri or fm
    scalar disMax = 0;
    label ncond = 0;

    //scalar disX = (maxInX - minInX)/ri_;
    //if(disX > disMax)
    //{
        //disMax = disX;
        //ncond = 0;
    //}

    //scalar disY = (maxInY - minInY)/ri_;
    //if(disY > disMax)
    //{
        //disMax = disY;
        //ncond = 1;
    //}

    //scalar disZ = (maxInZ - minInZ)/ri_;
    //if(disZ > disMax)
    //{
        //disMax = disZ;
        //ncond = 2;
    //}

    forAll(Xii_,i)
    {
        scalar disXiR = mag(maxInXiR[i] - minInXiR[i])/Xii_[i];
        if(disXiR > disMax)
        {
            disMax = disXiR;
            ncond = 3+i;
            //Info << "Updated ncond" << endl;
        }
    }


    lessArg comp(ncond);
    std::sort
    (
        iterL,
        iterU,
        [&](label& A, label& B) -> bool
        {
            return comp(particleList[A],particleList[B]);
        }
    );

    //- Recursive function calls for lower and upper branches of the particle list
    kdTreeLikeSearch(particleList,l,m,pInd,L,U);

    kdTreeLikeSearch(particleList,m+1,u,pInd,L,U);
};
