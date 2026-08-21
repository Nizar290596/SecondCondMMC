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
    Foam::ransMixParticleModel

Description
    Particle mixing model for RANS turbulence 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef ransMixParticleModel_H
#define ransMixParticleModel_H

#include "CloudMixingModel.H"
#include "baseParticleDataContainer.H"
#include "particlePtrList.H"
#include <vector>
#include <list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ransMixParticleModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ransMixParticleModel
:
    public CloudMixingModel<CloudType>
{
    typedef typename CloudType::particleType particleType;

    protected:
        
        // Struct to store the eulerian field properties for each particle
        class eulerianFieldData
        : public baseParticleDataContainer
        {            
            public:
                eulerianFieldData()
                : baseParticleDataContainer(2,3,1) {};
                
                //- Construct from input stream
                eulerianFieldData(Istream& is)
                : baseParticleDataContainer(is) 
                {};
                    
            // Non-Const Access
                //- random number
                inline scalar& Rand() {return scalarValues_[0];}
                
                //- tau turb
                inline scalar& tauTurb() {return scalarValues_[1];}
                
                //- Effective diffusivity
                inline scalar& DEff() {return scalarValues_[2];}
                
                //- Index of the particleList_
                inline label& particleIndex() {return labelValues_[0];}
        
                inline label& processorIndex() {return labelValues_[1];}
                
                inline scalarField& XiR() {return scalarFields_[0];}

            // Const Access
                //- random number
                inline const scalar& Rand() const {return scalarValues_[0];}
                inline const scalar& tauTurb() const {return scalarValues_[1];}
                inline const scalar& DEff() const {return scalarValues_[2];}
                inline const label& particleIndex() const {return labelValues_[0];}
                inline const label& processorIndex() const {return labelValues_[1];}
                inline const scalarField& XiR() const {return scalarFields_[0];}
        };
        
        //- Mixing time scale
        const mmcVarSet& XiR_;

        const wordList XiRNames_;
        
        //- Indexes in ransMixParticleModel
        const label numXiR_;

        //- Turbulent time scale 
        const volScalarField& tauTurb_;

        //- List of all particles date fields for each super cell 
        List<particlePtrList<particleType>> particleList_;
        
        //- List of all eulerian date fields for each super cell 
        List<DynamicList<eulerianFieldData>> eulerianFieldsList_;

        
    protected:
        
        //- Build particle list
        void buildParticleList();

        virtual void SmixList
        (
            particlePtrList<particleType>& particleList,
            DynamicList<eulerianFieldData>& eulFields
        ) = 0;
 
    public:
    // Constructors

        //- Construct from dictionary
        ransMixParticleModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type,
            const mmcVarSet& Xi
        );

        //- Construct copy
        ransMixParticleModel 
        (
            const ransMixParticleModel<CloudType>& cm
        );
        

    // Destructor
        virtual ~ransMixParticleModel();


    // Member functions

        //- Apply mix
        virtual void Smix();
        
    // Inline
        
        inline label numXiR() const
        {
            return numXiR_;
        };
        
        inline const wordList& XiRNames() const
        {
            return XiRNames_;
        };

        inline const mmcVarSet& XiR() const
        {
            return XiR_;
        };

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ransMixParticleModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
