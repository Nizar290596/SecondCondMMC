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
    Foam::kdTreeLikeMatching
 
Description
    Generate particle matching by constructing a k-d tree with box size 2. 
    The particles in the leaf box are then selected as pairs.


SourceFiles
    kdTreeLikeMatching.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/



#ifndef kdTreeLikeMatching_H
#define kdTreeLikeMatching_H

#include "dictionary.H"
#include "Tuple2.H"
#include "List.H"
#include "DynamicList.H"
#include "particleMatchingAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class kdTreeLikeMatching Declaration
\*---------------------------------------------------------------------------*/
template<class particleType>
class kdTreeLikeMatching
:
    public particleMatchingAlgorithm<particleType>
    //public particleMatchingAlgorithm<eulerianFieldData>
{
    private:

        //- Normalisation factors
        const scalar ri_;           // for physical coordinates
        const List<scalar> Xii_;    // for reference variable

        //- Search for particle pairs by using a recursive call to a 
        //- kdTree construction with box size 2
        void kdTreeLikeSearch
        (
            const DynamicList<particleType>& particleList,
            //const DynamicList<eulerianFieldData>& particleList,
            label l,
            label u,
            std::vector<label>& pInd,
            std::vector<label>& L,
            std::vector<label>& U
        ) const;

        // Required for kd-Tree like construction of finding particle pairs
        class lessArg
        {
            label ii_=-1;
            
            public:
        
                lessArg(label ii) {ii_=ii;} 

                inline bool operator()  
                (       
                    const particleType& one,
                    const particleType& two
                    //const eulerianFieldData& one,
                    //const eulerianFieldData& two
                ) const
                {
                    if (ii_ < 3)
                    {
			//Info << "Sorting in physical space " << endl;
                        return one.position()[ii_] < two.position()[ii_]; 
                    }
                    else 
                    {
			//Info << "Sorting in reference space " << endl;
                        label k = ii_-3;
                        return one.XiR()[k]< two.XiR()[k];
                    }

                    // Return a default to avoid warnings of the compiler
                    FatalError << "Comparison operator did not work properly"<<exit(FatalError);
                    return false;
                };
        };

    public:

        //- Construct with normalisation factor
        kdTreeLikeMatching
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
}

#include "kdTreeLikeMatching.C"

#endif

