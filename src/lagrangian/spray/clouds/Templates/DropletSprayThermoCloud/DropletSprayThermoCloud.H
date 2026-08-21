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
    Foam::DropletSprayThermoCloud

Description
    Templated base class for a thermo cloud with droplets coupling to mmcFoam
    for mass and energy exchange terms.

    <CloudName>Properties File:

    ... // Default Kinematic Cloud Options

    subModels
    {
        // Coupling model of droplets to stochastic particles
        dropletToMMCModel nearestNeighbor;
        nearestNeighborCoeffs {}

        // submodel for the additional fluid properties, such as
        // the binary diffusion coefficient
        liquidProperties UDF;
        UDFCoeffs {}

        dropletComponent    C2H5OH;     // Name of the droplet component/species
    
    }


SourceFiles
    DropletSprayThermoCloudI.H
    DropletSprayThermoCloud.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef DropletSprayThermoCloud_H
#define DropletSprayThermoCloud_H

#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "thermoCloud.H"
#include "baseParticleThermoMixture.H"

// Sub-models
#include "DropletToMMCModel.H"
#include "LiquidPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Create an enum type for the integration scheme
namespace Foam
{
    namespace DropletSprayThermo
    {
        enum class integrationScheme
        {
            EulerExplicit,
            CrankNicolson
        };
    }
}


namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DropletSprayThermoCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class mmcCloudType>
class DropletSprayThermoCloud
:
    public CloudType,
    public thermoCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef DropletSprayThermoCloud<CloudType, mmcCloudType> thermoCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<DropletSprayThermoCloud<CloudType,mmcCloudType>> cloudCopyPtr_;

        //- Access to the thermo model 
        //  Note: This allows access to the functions without modifying the 
        //        state of the mixture. Otherwise the carrier() function 
        //        of SLGThermo can be used. 
        autoPtr<baseParticleThermoMixture> thermoModels_;

        DropletSprayThermo::integrationScheme integrationScheme_;

    // Private member functions

        //- No copy construct
        DropletSprayThermoCloud(const DropletSprayThermoCloud&) = delete;

        //- No copy assignment
        void operator=(const DropletSprayThermoCloud&) = delete;


protected:


    // Proctected Member Variables

        refPtr<mmcCloudType> mmcCloud_;

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;
 
        //- Name of the droplet species -- only single component
        word speciesName_;
 
        // References to the carrier gas fields
 
            //- SLG thermodynamics package
            const SLGThermo& thermo_;

            //- Temperature [K]
            const volScalarField& T_;

            //- Pressure [Pa]
            const volScalarField& p_;

            //- Mass transfer field for the droplet species
            //  The droplet can only have a single component! 
            autoPtr<volScalarField::Internal> rhoTrans_;

        // References to the cloud sub-models

            //- Coupling model from the droplets to the stochastic particles
            autoPtr
            <
                DropletToMMCModel
                <DropletSprayThermoCloud<CloudType,mmcCloudType>, mmcCloudType>
            >
                dropletToMMCModel_;

            autoPtr<LiquidPropertiesModel<thermoCloudType>> 
                liquidProperties_;

    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();

        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(DropletSprayThermoCloud<CloudType,mmcCloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        DropletSprayThermoCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        DropletSprayThermoCloud
        (
            DropletSprayThermoCloud<CloudType,mmcCloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        DropletSprayThermoCloud
        (
            const fvMesh& mesh,
            const word& name,
            const DropletSprayThermoCloud<CloudType,mmcCloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new DropletSprayThermoCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new DropletSprayThermoCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~DropletSprayThermoCloud() = default;

    // Member Functions
 
        // Access
 
            //- Return a reference to the cloud copy
            inline const thermoCloudType& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;
 
            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();
 
            //- Return const access to thermo package
            inline const SLGThermo& thermo() const;
 
            //- Access the thermo models, e.g., transport, equation of state,...
            const baseParticleThermoMixture& thermoModels() const;

            //- Additional fluid properties, e.g., sat pressure, 
            //  diffusion coefficient, latent heat of vaporization etc. 
            const LiquidPropertiesModel<thermoCloudType>& 
                liquidProperties() const;

            //- Return const access to the carrier temperature field
            inline const volScalarField& T() const;
 
            //- Return const access to the carrier pressure field
            inline const volScalarField& p() const;

            //- Return the species name 
            inline const word& dropletSpeciesName() const;

        // Sources

            //- Return const access to mass source fields
            inline const volScalarField::Internal&
                rhoTrans() const;

            //- Return reference to mass source fields
            inline volScalarField::Internal&
                rhoTrans();

            //- Return tmp total mass source for carrier phase
            //  - fully explicit
            inline tmp<volScalarField::Internal> Srho() const;

            //- Return total mass source term [kg/m3/s]
            inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;

        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources
            (
                const DropletSprayThermoCloud<CloudType,mmcCloudType>& cloudOldTime
            );

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve(const typename parcelType::trackingData& td);

            //- Evolve the cloud
            void evolve();

        // Integration Scheme Selection

            //- Return the enum type of the integration scheme
            inline const DropletSprayThermo::integrationScheme& 
            integrationScheme() const;
            


        // MMC Reference

            //- Set the mmc cloud
            void setMMCCloud(mmcCloudType& cloud);

            //- Const access the mmc cloud
            const mmcCloudType& mmcCloud() const;

            //- Access the mmc cloud
            mmcCloudType& mmcCloud();

            //- Access to the coupling model
            const DropletToMMCModel
            <
                DropletSprayThermoCloud<CloudType,mmcCloudType>,
                mmcCloudType
            >&
            dropletToMMCCoupling() const;


        // Radiation - Currently returning blank fields
 
            //- Return tmp equivalent particulate emission
            virtual tmp<volScalarField> Ep() const;

            //- Return tmp equivalent particulate absorption
            virtual tmp<volScalarField> ap() const;

            //- Return tmp equivalent particulate scattering factor
            virtual tmp<volScalarField> sigmap() const;


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DropletSprayThermoCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DropletSprayThermoCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
