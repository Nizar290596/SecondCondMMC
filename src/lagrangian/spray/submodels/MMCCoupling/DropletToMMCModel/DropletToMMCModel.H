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
    Foam::DropletToMMCModel

Description
    Base class for coupling droplets of the spray to the mmc cloud

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        dropletToMMCCoupling    equal;
    }
    \endverbatim

SourceFiles
    DropletToMMCModel.C
    DropletToMMCModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DropletToMMCModel_H
#define DropletToMMCModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DropletToMMCModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class mmcCloudType>
class DropletToMMCModel
:
    public CloudSubModelBase<CloudType>
{
protected:

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

public:

    //- Runtime type information
    TypeName("DropletToMMCModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DropletToMMCModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Generated Methods

        //- No copy assignment
        void operator=(const DropletToMMCModel<CloudType, mmcCloudType>&) = delete;


    // Constructors

        //- Construct null from owner
        DropletToMMCModel(CloudType& owner);

        //- Construct from dictionary
        DropletToMMCModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Copy construct
        DropletToMMCModel(const DropletToMMCModel<CloudType, mmcCloudType>& htm);

        //- Construct and return a clone
        virtual autoPtr<DropletToMMCModel<CloudType, mmcCloudType>> clone() const = 0;


    //- Destructor
    virtual ~DropletToMMCModel() = default;


    //- Selector
    static autoPtr<DropletToMMCModel<CloudType, mmcCloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Evaluation

            //- Couple the mass and enthalpy sources to the mmc stochastic
            //  particle
            //  This function is not const, as the stochastic particles in the
            //  mmcCloud are modified
            //  Evaporation/Condensation source term [kg/s] -- mDot
            //  Enthalpy exchange term d(m h)/dt [J/s]      -- hDot
            virtual void coupleSources
            (
                const parcelType& p,    // Particle
                const typename parcelType::trackingData& td,
                mmcCloudType& mmcCloud,
                const scalar& mDot,
                const scalar& hDot,
                const label jFuel,
                const scalar& dt
            ) const = 0;


            //- Collect the information from the stochastic particles 
            //  to get the gas phase conditions
            virtual void getGasPhaseConditions
            (
                const parcelType& p,    // Particle
                const typename parcelType::trackingData& td,
                mmcCloudType& mmcCloud,
                scalar& TInf,
                scalar& pInf,
                List<scalar>& YInf
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDropletToMMCModel(CloudType, mmcCloudType)                         \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
                                                                               \
    typedef DropletToMMCModel<thermoCloudType,mmcCloudType>                    \
        dropletToMMCModelType;                                                 \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        dropletToMMCModelType,                                                 \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            dropletToMMCModelType,                                             \
            dictionary                                                         \
        );                                                                     \
    }


#define makeDropletToMMCModelType(SS, CloudType, mmcCloudType)                 \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    typedef Foam::SS<thermoCloudType,mmcCloudType>                             \
        dropletToMMCModelSSType;                                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
        (dropletToMMCModelSSType, 0);                                          \
                                                                               \
    Foam::DropletToMMCModel<thermoCloudType,mmcCloudType>::                    \
        adddictionaryConstructorToTable                                        \
            <dropletToMMCModelSSType>                                          \
            add##SS##CloudType##thermoCloudType##mmcCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DropletToMMCModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
