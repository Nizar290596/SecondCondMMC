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
    Foam::OFLiquidProperties

Description
    Get liquid properties form OpenFOAM databank.
    Uses liquidProperties class 

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        saturationPressure    OFLiquidProperties;

        OFCoeffs
        {
            fluidName   "Ethanol";
        }
    }
    \endverbatim

SourceFiles
    OFLiquidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef OFLiquidProperties_H
#define OFLiquidProperties_H

#include "LiquidPropertiesModel.H"
#include "Function1.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class OFLiquidProperties Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class OFLiquidProperties
:
    public LiquidPropertiesModel<CloudType>
{

private:

    // Private Variables

        autoPtr<liquidProperties> liqProp_;

public:

    //- Runtime type information
    TypeName("OFLiquidProperties");


    // Generated Methods

        //- No copy assignment
        void operator=(const OFLiquidProperties&) = delete;


    // Constructors

        //- Construct from dictionary
        OFLiquidProperties(const dictionary&, CloudType& owner);

        //- Copy construct
        OFLiquidProperties(const OFLiquidProperties<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<LiquidPropertiesModel<CloudType>> clone() const
        {
            return autoPtr<LiquidPropertiesModel<CloudType>>
            (
                new OFLiquidProperties<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~OFLiquidProperties() = default;


    // Member Functions

        //- Calculate the saturation pressure based on the given temperature
        virtual scalar pSat(const scalar T) const;

        //- Calculate the saturation temperature
        virtual scalar TSat(const scalar p) const;


        //- Calculate the diffusion coefficient
        virtual scalar Df
        (
            const scalar p,
            const scalar T,
            const scalar Wb
        ) const;

        //- Calculate the latent heat
        virtual scalar Lv(const scalar p, const scalar T) const;

        //- Calculate the heat capacity of the liquid
        virtual scalar CpL(const scalar p, const scalar T) const;

        //- Calculate the liquid density
        virtual scalar rhoL(const scalar p, const scalar T) const;

        //- Calculate the absolute enthalpy of the liquid
        virtual scalar Ha(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OFLiquidProperties.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
