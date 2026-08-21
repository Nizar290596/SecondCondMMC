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
    along with  mmcFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Test the kdTree search algorithm in comparison to a brute force approach
    on a real data set

    Times of construction and search are reported

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "kdTree.H"
#include "fvCFD.H"
#include "IFstream.H"
#include <chrono>
#include <random>
#include <queue>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace kdTreeTest
{
    Foam::DynamicList<Foam::List<Foam::label>> findPairs
    (
        const List<List<scalar>>& particles,
        const List<label>& dims,
        const List<scalar>& weigths
    );
    
    void kdTreeLikeSearch
    (
        const List<List<scalar>>& particleList,
        const List<scalar>& weigths,
        label l,
        label u,
        std::vector<label>& pInd,
        std::vector<label>& L,
        std::vector<label>& U    
    );
    
}

TEST_CASE("k-d Tree Pairing Test","[kdTree][mmcSupport]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();

    // DATA SETTINGS

        // Data for dimensions and weights taken from an example Shear Layer test
        // case
        labelList dims(4);
        dims[0] = 1;
        dims[1] = 2;
        dims[2] = 3;
        dims[3] = 13;

        scalarList wts(4);
        wts[0] = 1;
        wts[1] = 1;
        wts[2] = 1;
        wts[3] = 0.03; // 0.03;
        
    // LOAD INPUT DATA

        // First load the list of particles from the case
        Foam::IFstream ifs(args.path()/"particleDataSets/particleList.dat");
        const DynamicList<List<scalar>> particleListFull(ifs);
        
        // Search 20 nearest neighbours for each cell 
        // Read the query points from file
        Foam::IFstream ifsQueryPoints(args.path()/"particleDataSets/queryPointList.dat");
        const List<scalarList> queryPoints(ifsQueryPoints);
        const label nQueryPoints = queryPoints.size();

        // Create a particle List based on the given dims
        List<List<scalar>> particleList(particleListFull.size());
        forAll(particleList,i)
        {
            List<scalar> temp(4);
            for (label k=0; k < 4; k++)
                temp[k] = particleListFull[i][dims[k]];
            particleList[i] = temp;
        }


    // Terminal Output
    Info << "Finding Particle Pairs:"<<endl;


    // Find pairs 

        // Find the unique pairs
        
        auto t1 = std::chrono::high_resolution_clock::now();
        kdTree<List<scalar>,2> tree(particleList,wts,true);
        auto pairs = tree.leafParticles();
        auto t2 = std::chrono::high_resolution_clock::now();
        
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        Info << "\tFinding pairs took: " << ms_int.count() << "ms" << endl;
        
        Info << "\tOf nParticles "<<particleList.size()
             << "   #"<<pairs.size()<<" pairs have been found"<<endl;
             
        // Check that the number of particles in each bucket is not larger 
        // than 3
        label nPairsWithThreeEntries = 0;
        for (auto& pair : pairs)
        {
            REQUIRE(pair.size() <= 3);
            REQUIRE(pair.size() > 1);
            if (pair.size() == 3)
                nPairsWithThreeEntries++;
        }
        Info << "\tNumber of pairs with three entries: "<<nPairsWithThreeEntries<<endl;
        
        // ====================================================================
        // Compare to search algorithm as it is currently implemented
        Info << "Searching particle pairs with old k-d tree like search"<<endl;
        auto t3 = std::chrono::high_resolution_clock::now();
        auto pairs2 = kdTreeTest::findPairs
        (
            particleList,dims,wts
        );
        auto t4 = std::chrono::high_resolution_clock::now();
        
        auto ms_int_test = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
        Info << "\tFinding pairs old-kd-Tree took: " << ms_int_test.count() << "ms" << endl;
        Info << "\tOf nParticles "<<particleList.size()
             << "   #"<<pairs2.size()<<" pairs have been found"<<endl;

        nPairsWithThreeEntries = 0;
        for (auto& pair : pairs2)
        {
            REQUIRE(pair.size() <= 3);
            REQUIRE(pair.size() > 1);
            if (pair.size() == 3)
                nPairsWithThreeEntries++;
        }
        
        Info << "\tNumber of pairs with three entries: "<<nPairsWithThreeEntries<<endl;
}


Foam::DynamicList<Foam::List<Foam::label>> kdTreeTest::findPairs
(
    const List<List<scalar>>& particles,
    const List<label>& dims,
    const List<scalar>& weights
) 
{   
    DynamicList<List<label>> pairs;
    // Keeping track of indices for kdTreeLikeSearch
    std::vector<label> L;
    std::vector<label> U;
    L.reserve(particles.size());
    U.reserve(particles.size());
    
    // create an index list for the particle data
    std::vector<label> pInd(particles.size());
    std::iota(pInd.begin(),pInd.end(),0);
    
    kdTreeLikeSearch(particles,weights,1,particles.size(),pInd,L,U);

    // reserve space for list of pairs
    pairs.reserve(ceil(0.5*particles.size()));

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
    
    return pairs;
}


void kdTreeTest::kdTreeLikeSearch
(
    const List<List<scalar>>& particleList,
    const List<scalar>& weights,
    label l,
    label u,
    std::vector<label>& pInd,
    std::vector<label>& L,
    std::vector<label>& U    
)
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
    scalar maxInXiR = 0;
    
    scalar minInX = GREAT;
    scalar minInY = GREAT;
    scalar minInZ = GREAT;
    scalar minInXiR = GREAT;
    
    // Find minimum and maximum for each coordinate
    for (auto it = iterL; it != iterU; it++)
    {
        auto& pos = particleList[*it];
        maxInX = std::max(maxInX,pos[0]);
        maxInY = std::max(maxInY,pos[1]);
        maxInZ = std::max(maxInZ,pos[2]);
        
        minInX = std::min(minInX,pos[0]);
        minInY = std::min(minInY,pos[1]);
        minInZ = std::min(minInZ,pos[2]);
    
        maxInXiR = std::max(maxInXiR,pos[3]);
        minInXiR = std::min(minInXiR,pos[3]);
    }

    //- Scaled/stretched distances between Max and Min in each direction
    //- Default is random mixing, overwritten if mixing distances greater than ri or fm
    scalar disMax = 0;
    label ncond = 0;

    scalar disX = (maxInX - minInX)/weights[0];
    if(disX > disMax)
    {
        disMax = disX;
        ncond = 0;
    }

    scalar disY = (maxInY - minInY)/weights[1];
    if(disY > disMax)
    {
        disMax = disY;
        ncond = 1;
    }

    scalar disZ = (maxInZ - minInZ)/weights[2];
    if(disZ > disMax)
    {
        disMax = disZ;
        ncond = 2;
    }

    scalar disXi = (maxInXiR - minInXiR)/weights[3];
    if(disXi > disMax)
    {
        disMax = disXi;
        ncond = 3;
    }


    auto comp = [&ncond](const List<scalar>& l1, const List<scalar>& l2)
    {
        return l1[ncond] < l2[ncond];
    };
    
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
    kdTreeLikeSearch(particleList,weights,l,m,pInd,L,U);

    kdTreeLikeSearch(particleList,weights,m+1,u,pInd,L,U);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

