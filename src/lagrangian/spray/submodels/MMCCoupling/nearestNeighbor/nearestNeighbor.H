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
    Foam::nearestNeighbor

Description
    Equal model for coupling spray droplets to the mmc stochastic particles.
    This model distributes the source to all stochastic particles in a given 
    volume. 

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        dropletToMMCCoupling    nearestNeighbor;
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      heatTransferModel  | Type name: none          | word | yes  | -
    \endtable

SourceFiles
    nearestNeighbor.C

\*---------------------------------------------------------------------------*/

#ifndef nearestNeighbor_H
#define nearestNeighbor_H

#include "DropletToMMCModel.H"
#include "kdTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nearestNeighbor Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class mmcCloudType>
class nearestNeighbor
:
    public DropletToMMCModel<CloudType, mmcCloudType>
{
protected:

    // Typedefs for the mmc stochastic particles
    typedef typename mmcCloudType::particleType mmcParticleType;
    typedef typename IDLList<mmcParticleType>::iterator mmcParticleIteratorType;

    //- Convenience typedef to the cloud's parcel type
    typedef typename CloudType::parcelType parcelType;

public:
    
    // Private struct required for finding the closest stochastic particles 
    // with the k-d tree tool
    class  container
    : public List<scalar>
    {
        public:

            // Iterator pointing to the particle 
            mmcParticleType* iter_;

        // Constructor

            //- Default
            container() = default;

            //- With entries
            container
            (
                bool lambda_x, bool lambda_f, bool lambda_T,
                scalar x, scalar y, scalar z, scalar f, scalar T, 
                mmcParticleType* iter
            )
            {
                label size = label(lambda_x)*3+label(lambda_T)+label(lambda_f);

                this->resize(size);
                label ind = 0;
                if (lambda_x)
                {
                    this->operator[](ind++) = x;
                    this->operator[](ind++) = y;
                    this->operator[](ind++) = z;
                }
                if (lambda_f)
                    this->operator[](ind++) = f;
                if (lambda_T)
                    this->operator[](ind++) = T;
                
                #ifdef FULLDEBUG
                    if (size != ind)
                        FatalError 
                            << "Size does not match in container of droplet to mmc coupling"
                            << exit(FatalError);
                #endif

                iter_ = iter;
            }
    };

private:

    // Private Member Variables

        //- normalization parameters
        const bool lambda_x;    // Switches on using the coordinates for distance
        const bool lambda_f;
        const bool lambda_T;

        const scalar ri_;
        const scalar fm_;
        const scalar Tm_;

        const bool useSPProps_;

        //- List of stochastic partices the kdTree refers to
        mutable DynamicList<container> particleList_;

        //- k-d tree to find the closest stochastic particle
        mutable autoPtr<kdTree<container>> kdTree_;

        //- time index of the current time step
        //  Required to determine if the k-d tree needs to be regenerated
        mutable label timeIndex_{-1};



    // Private Member Functions

        //- Generate the stochastic particle list for k-d tree generation
        void buildParticleList
        (
            mmcCloudType& mmcCloud
        ) const;

        //- Generate the k-d tree if necessary
        void genKdTree(mmcCloudType& mmcCloud) const;

        mmcParticleType* locateClosestStochasticParticle
        (
            const parcelType& p,    // Particle
            const typename parcelType::trackingData& td,
            mmcCloudType& mmcCloud
        ) const;



public:

    //- Runtime type information
    TypeName("nearestNeighbor");


    // Generated Methods

        //- No copy assignment
        void operator=(const nearestNeighbor&) = delete;


    // Constructors

        //- Construct from dictionary
        nearestNeighbor(const dictionary&, CloudType& owner);

        //- Copy construct
        nearestNeighbor(const nearestNeighbor<CloudType, mmcCloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<DropletToMMCModel<CloudType, mmcCloudType>> clone() const
        {
            return autoPtr<DropletToMMCModel<CloudType, mmcCloudType>>
            (
                new nearestNeighbor<CloudType, mmcCloudType>(*this)
            );
        }


    //- Destructor
    virtual ~nearestNeighbor() = default;


    // Member Functions

        //- Couple the mass and enthalpy sources to the mmc stochastic
        //  particle
        //  Requires the position of the droplet        -- pos
        //  Reference variables of the droplet          -- refVars
        //  Evaporation/Condensation source term [kg/s] -- mDot
        //  Enthalpy exchange term d(m h)/dt [J/s]      -- hDot
        void coupleSources
        (
            const parcelType& p,    // Particle
            const typename parcelType::trackingData& td,
            mmcCloudType& mmcCloud,
            const scalar& mDot,
            const scalar& hDot,
            const label jFuel,
            const scalar& dt
        ) const;

        //- Collect the information from the stochastic particles 
        //  to get the gas phase conditions
        void getGasPhaseConditions
        (
            const parcelType& p,    // Particle
            const typename parcelType::trackingData& td,
            mmcCloudType& mmcCloud,
            scalar& TInf,
            scalar& pInf,
            List<scalar>& YInf
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nearestNeighbor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
