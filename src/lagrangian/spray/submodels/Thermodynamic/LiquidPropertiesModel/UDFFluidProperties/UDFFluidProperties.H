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
    Foam::UDFFluidProperties

Description
    User defined function for the saturation pressure. 
    Utilizes the Function1<scalar> concept of OpenFOAM.

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        saturationPressure    UDF;

        UDFCoeffs
        {
            saturationPressure 
            {
                type    constant 10;     // can also use coded
            }

            saturationTemperature
            {
                type    constant 10;     // can also use coded
            }

            binaryDiffusionCoefficient 1E-5;

            latentHeat
            {
                type        table ();
            }

            liquidHeatCapacity
            {
                // Note: for coded only the temperature input is considered
                type        coded;
                code
                #{
                    // Note: identifier x is the input variable from the function
                    scalar a = 1.0;
                    return a;
                #};
            }

            liquidAbsoluteEnthalpy
            {
                type        table ();   // You can also use coded or constant,
                                        // or any other Function1<> object type
            }

            // Liquid density 
            rhoL
            {
                type        table ();   // You can also use coded or constant,
                                        // or any other Function1<> object type
            }
        }
    }
    \endverbatim

SourceFiles
    UDFFluidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef UDFFluidProperties_H
#define UDFFluidProperties_H

#include "LiquidPropertiesModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class UDFFluidProperties Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class UDFFluidProperties
:
    public LiquidPropertiesModel<CloudType>
{

private:

    // Private Variables

        //- General purpose Function1 object for the saturation pressure
        autoPtr<Function1<scalar>> saturationPressure_;

        //- General purpose Function1 object for the saturation temperature
        autoPtr<Function1<scalar>> saturationTemperature_;

        //- General purpose Function1 object for the latent heat
        //  returns the latent heat of vaporization as [J/kg]
        autoPtr<Function1<scalar>> latentHeat_;

        //- General purpose Function1 object for the heat capacity
        //  returns the heat capacity 
        autoPtr<Function1<scalar>> heatCapacityLiquid_;

        //- General purpose Function1 object for the absolute enthalpy.
        autoPtr<Function1<scalar>> absoluteEnthalpy_;

        //- General purpose Function1 object for the liquid density
        autoPtr<Function1<scalar>> rhoL_;

        //- Binary diffusio coefficient
        //  Uses the Fuller-Schettler-Giddings equation
        //  and combines all constant values into one coefficient
        //  D_A-air = Coeff T^{1.75}/p;
        const scalar binaryDiffusionCoefficient_;

public:

    //- Runtime type information
    TypeName("UDF");


    // Generated Methods

        //- No copy assignment
        void operator=(const UDFFluidProperties&) = delete;


    // Constructors

        //- Construct from dictionary
        UDFFluidProperties(const dictionary&, CloudType& owner);

        //- Copy construct
        UDFFluidProperties(const UDFFluidProperties<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<LiquidPropertiesModel<CloudType>> clone() const
        {
            return autoPtr<LiquidPropertiesModel<CloudType>>
            (
                new UDFFluidProperties<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~UDFFluidProperties() = default;


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
    #include "UDFFluidProperties.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
