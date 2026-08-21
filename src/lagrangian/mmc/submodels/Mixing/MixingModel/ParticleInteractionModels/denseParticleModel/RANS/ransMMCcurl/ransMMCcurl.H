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
    Foam::ransMMCcurl

Description
    RANS-MMC mixing based on Curl's model
     
    Note: For particle pair searching only the distance in the reference 
          variable space is used!

\*---------------------------------------------------------------------------*/

#ifndef ransMMCcurl_H
#define ransMMCcurl_H

#include "ransMixParticleModel.H"
#include <list>
#include "kdTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ransMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ransMMCcurl
:
    public ransMixParticleModel<CloudType>
{
    typedef typename CloudType::particleType particleType;
    typedef typename ransMixParticleModel<CloudType>::eulerianFieldData 
        eulerianFieldData;
   
    // Private data
    
        //- Time scale model coefficients
        const scalar C_phi_; // Unconditional scalar dissipation constant
        const scalar rt_; // Target correlation coefficient
        const scalar bo_; 

        const List<scalar> Xi_i_; // Normalization factor

        const mmcVarSet& XiR_;
    
        // Required for kd-Tree like construction of finding particle pairs
        class lessArg
        {
            label ii_=-1;
            
            public:
        
                lessArg(label ii) {ii_=ii;} 

                inline bool operator()  
                (       
                    const eulerianFieldData& one,
                    const eulerianFieldData& two
                ) const
                {
                    if (ii_ < 3)
                    {
                        return one.position()[ii_] < two.position()[ii_]; 
                    }
                    else 
                    {
                        label k = ii_-3;
                        return one.XiR()[k]< two.XiR()[k];
                    }

                    // Return a default to avoid warnings of the compiler
                    FatalError << "Comparison operator did not work properly"<<exit(FatalError);
                    return false;
                };
        };
    



    protected:
    
    
        //- Mixing particle pairs    
        void mixpair
        (
            particleType& p, 
            const eulerianFieldData& pEulField,
            particleType& q,
            const eulerianFieldData& qEulField,
            scalar& deltaT
        );
        
        void mixScalarProperties
        (
            particleType& p,
            particleType& q,
            scalar mixExtent
        );        
        
        virtual void SmixList
        (
            particlePtrList<particleType>& particleList,
            DynamicList<eulerianFieldData>& eulFields
        );
    
        virtual void MixXi
        (
            particlePtrList<particleType>& particleList,
            DynamicList<eulerianFieldData>& eulFields
        );
        
        //- Search for particle pairs by using a recursive call to a 
        //- kdTree construction with box size 2
        //  NOTE: For the RANS model only the distance in reference space is
        //        used
        void kdTreeLikeSearch
        (
            label l,
            label u,
            DynamicList<eulerianFieldData>& particleList,
            std::vector<label>& L,
            std::vector<label>& U
        );

        //- Read the normalisation factors for the reference variable 
        //- from the dictionary
        List<scalar> getXiNormalisation();
    
    public:

        //- Runtime type information
        TypeName("ransMMCcurl");

    // Constructors

        //- Construct from dictionary
        ransMMCcurl(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        ransMMCcurl(const ransMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new ransMMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~ransMMCcurl();
        
    // Member Function
        
        virtual const scalarField XiR0(label patch, label patchFace);
        virtual const scalarField XiR0(label celli);      

        //- return mmcVarSet
        virtual const mmcVarSet& XiR() const;
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ransMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
