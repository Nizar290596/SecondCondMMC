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
    Foam::denseMixParticleModel

Description
    Mixing model for dense LES approach
    Super mesh and LES mesh have to be equal!

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef denseMixParticleModel_H
#define denseMixParticleModel_H

#include "CloudMixingModel.H"
#include <vector>
#include <list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class denseMixParticleModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class denseMixParticleModel
:
    public CloudMixingModel<CloudType>
{
    typedef typename CloudType::particleType particleType;
    
    public:
        // Struct to store the eulerian field properties for each particle
        class eulerianFieldParticleValues
        {
            private:
                //- Fixed List to store the eulerian field values for each particle
                FixedList<scalar,3> values_;
            
            public:
                eulerianFieldParticleValues() = default;
                    
                // Non-Const Access
                    //- random number
                    inline scalar& Rand() {return values_[0];}
                    
                    //- Squared filter width
                    inline scalar& sqrFilter() {return values_[1];}
                    
                    //- Effective diffusivity
                    inline scalar& DEff() {return values_[2];}

                // Const Access
                    //- random number
                    inline const scalar& Rand() const {return values_[0];}
                    
                    //- Effective diffusivity 
                    inline const scalar& sqrFilter() const {return values_[1];}
                    
                    //- Effective diffusivity 
                    inline const scalar& DEff() const {return values_[2];}
        };
    
    protected:
    
        //- sqr of filter width and turbulence field parameters
        const volScalarField& DEff_;

        //- List of particles in each LES cell
        List<DynamicList<particleType*>> particleList_;

        //- List of all eulerian date fields for each particle for each LES cell
        List<DynamicList<eulerianFieldParticleValues>> eulerianFieldsList_;

        //- Mix particle for each list of particles in each cell
        virtual void SmixList
        (
            DynamicList<particleType*>& pList,
            DynamicList<eulerianFieldParticleValues>& eulFields
        ) = 0;

    public:
    // Constructors

        //- Construct from dictionary
        denseMixParticleModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        denseMixParticleModel 
        (
            const denseMixParticleModel<CloudType>& cm
        );
        

    // Destructor
        virtual ~denseMixParticleModel();


    // Member functions
    
        //- Build particle list
        void buildParticleList();

        //- Apply mix
        virtual void Smix() override;
        
           
    // Inline
        
        inline label numXiR() const
        {
            return 0;
        };
        
        inline const wordList& XiRNames() const
        {
            return wordList::null();
        };

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "denseMixParticleModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
