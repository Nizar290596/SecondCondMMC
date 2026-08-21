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
    Foam::Curl

Description
    Mixing model for dense LES approach
    Mixing time based on Curl's model

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef Curl_H
#define Curl_H

#include "denseMixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Curl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Curl
:
    public denseMixParticleModel<CloudType>
{   
        typedef typename CloudType::particleType particleType;
        typedef typename denseMixParticleModel<CloudType>::eulerianFieldParticleValues
            eulerianFieldParticleValues;

    // Private data
    
        //- Time scale model coefficients
        const scalar CL_;
        const scalar CE_;

        //- Information about the reference variable 
        const mmcVarSet& XiR_;     
        
    protected:
        
        virtual void SmixList
        (
            DynamicList<particleType*>& particleList,
            DynamicList<eulerianFieldParticleValues>& eulFields
        );
        
        void mixpair
        (
            particleType& p,
            const eulerianFieldParticleValues& pEulField,
            particleType& q,
            const eulerianFieldParticleValues& qEulField,
            scalar& deltaT
        );
        
        void mixScalarProperties
        (
            particleType& p,
            particleType& q,
            scalar mixExtent
        );

    public:

        //- Runtime type information
        TypeName("Curl");

    // Constructors

        //- Construct from dictionary
        Curl(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        Curl(const Curl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new Curl<CloudType>(*this)
            );
        }

    // Destructor
        virtual ~Curl();

    // Member Functions
        
        //- Return scalar field with reference variables
        virtual inline const scalarField XiR0(label patch, label patchFace);
        virtual inline const scalarField XiR0(label celli);

        //- return mmcVarSet
        virtual const mmcVarSet& XiR() const;

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Curl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
