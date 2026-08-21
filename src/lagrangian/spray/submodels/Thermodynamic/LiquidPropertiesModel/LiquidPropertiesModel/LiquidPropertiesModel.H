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
    Foam::LiquidPropertiesModel

Description
    Base class for getting the liquid properties, e.g. density, 
    absolute enthalpy, but also saturation pressure and temperature 

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        additionalFluidProperties    UDF;
    }
    \endverbatim

SourceFiles
    LiquidPropertiesModel.C
    LiquidPropertiesModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef LiquidPropertiesModel_H
#define LiquidPropertiesModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class LiquidPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LiquidPropertiesModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("LiquidPropertiesModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        LiquidPropertiesModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Generated Methods

        //- No copy assignment
        void operator=(const LiquidPropertiesModel<CloudType>&) = delete;


    // Constructors

        //- Construct null from owner
        LiquidPropertiesModel(CloudType& owner);

        //- Construct from dictionary
        LiquidPropertiesModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Copy construct
        LiquidPropertiesModel(const LiquidPropertiesModel<CloudType>& satModel);

        //- Construct and return a clone
        virtual autoPtr<LiquidPropertiesModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~LiquidPropertiesModel() = default;


    //- Selector
    static autoPtr<LiquidPropertiesModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Calculate the saturation pressure based on the given temperature
        virtual scalar pSat(const scalar T) const = 0;

        //- Calculate the saturation temperature
        virtual scalar TSat(const scalar p) const = 0;

        //- Calculate the diffusion coefficient for given weight of the surrounding
        //  gas
        virtual scalar Df
        (
            const scalar p,
            const scalar T,
            const scalar Wb
        ) const = 0;

        //- Calculate the latent heat
        virtual scalar Lv(const scalar p, const scalar T) const = 0;

        //- Calculate the heat capacity of the liquid
        virtual scalar CpL(const scalar p, const scalar T) const = 0;

        //- Calculate the liquid density
        virtual scalar rhoL(const scalar p, const scalar T) const = 0;

        //- Calculate the absolute enthalpy of the liquid
        virtual scalar Ha(const scalar p, const scalar T) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLiquidPropertiesModel(CloudType)                          \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::LiquidPropertiesModel<thermoCloudType>,                      \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            LiquidPropertiesModel<thermoCloudType>,                        \
            dictionary                                                         \
        );                                                                     \
    }


#define makeLiquidPropertiesModelType(SS, CloudType)                  \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);         \
                                                                               \
    Foam::LiquidPropertiesModel<thermoCloudType>::                         \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType>>             \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LiquidPropertiesModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
