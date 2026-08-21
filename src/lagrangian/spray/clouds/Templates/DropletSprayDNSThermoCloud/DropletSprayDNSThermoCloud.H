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
    Foam::DropletSprayDNSThermoCloud

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

        sourceDistribution     PSI;     // Default is particle-source-in-cell model (PSI)
                                        // (Optional) activate distributed

        // If distSourceModel is active the sourceDistrbutionCoeff need 
        // to provide the averaging volume
        distributedCoeffs 
        {
            // Type of the averaging volume, either box or sphere        
            type       box;

            // Potentially add additional entries required for the averaging
            // volume
        }
    
    }


SourceFiles
    DropletSprayDNSThermoCloudI.H
    DropletSprayDNSThermoCloud.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef DropletSprayDNSThermoCloud_H
#define DropletSprayDNSThermoCloud_H

#include "KinematicCloud.H"
#include "SLGThermo.H"
#include "thermoCloud.H"
#include "baseParticleThermoMixture.H"
#include "filterDNSField.H"

// Sub-models
#include "DropletToMMCModel.H"
#include "LiquidPropertiesModel.H"
#include "distributeSource/distributeSource.H"
#include "movingAverage/movingAverage.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create an enum type for the integration scheme
namespace Foam
{
    namespace DropletSprayDNSThermo
    {
        enum class integrationScheme
        {
            EulerExplicit,
            CrankNicolson
        };

        enum class distributeSourceModel
        {
            PSI,
            distributed
        };
    }
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DropletSprayDNSThermoCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class mmcCloudType>
class DropletSprayDNSThermoCloud
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
        typedef DropletSprayDNSThermoCloud<CloudType,mmcCloudType> thermoCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<DropletSprayDNSThermoCloud<CloudType,mmcCloudType>> cloudCopyPtr_;

        //- Access to the thermo model 
        //  Note: This allows access to the functions without modifying the 
        //        state of the mixture. Otherwise the carrier() function 
        //        of SLGThermo can be used. 
        autoPtr<baseParticleThermoMixture> thermoModels_;

        //- Distribute Source function
        autoPtr<distributeSource> distSource_;
        autoPtr<movingAverage> movingAverage_;

        //- Type of the distribute source model
        DropletSprayDNSThermo::distributeSourceModel distSourceModelType_;

        DropletSprayDNSThermo::integrationScheme integrationScheme_;

    // Private member functions

        //- No copy construct
        DropletSprayDNSThermoCloud(const DropletSprayDNSThermoCloud&) = delete;

        //- No copy assignment
        void operator=(const DropletSprayDNSThermoCloud&) = delete;


protected:


    // Proctected Member Variables

        refPtr<mmcCloudType> mmcCloud_;

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;
 
        //- Name of the droplet species -- only single component
        word speciesName_;

        //- Mixing rule fraction for the reference conditions
        //  e.g., Tref = Ts + A*(TInf - Ts)
        //  with Ts as surface temperature and TInf as infinity gas temperature
        const scalar refPropertiesFactor_;
 
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

            //- Enthalpy field for the droplet species
            //  The droplet can only have a single component! 
            autoPtr<volScalarField::Internal> hsTrans_;

            //- Enthalpy Coeff. field for the droplet species
            //  The droplet can only have a single component! 
            autoPtr<volScalarField::Internal> hsCoeff_;

        // References to the cloud sub-models

            //- Coupling model from the droplets to the stochastic particles
            autoPtr
            <
                DropletToMMCModel
                <DropletSprayDNSThermoCloud<CloudType,mmcCloudType>, mmcCloudType>
            >
                dropletToMMCModel_;

            //- Coupling model from the droplets to the stochastic particles
            autoPtr<LiquidPropertiesModel<thermoCloudType>>
                liquidProperties_;

    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();

            //- Set the distribute source model
            void setDistributeSource
            (
                const fvMesh& mesh,
                const dictionary& dict
            );

        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(DropletSprayDNSThermoCloud<CloudType, mmcCloudType>& c);

public:

    // Constructors

        //- Construct given carrier gas fields
        DropletSprayDNSThermoCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        DropletSprayDNSThermoCloud
        (
            DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        DropletSprayDNSThermoCloud
        (
            const fvMesh& mesh,
            const word& name,
            const DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new DropletSprayDNSThermoCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new DropletSprayDNSThermoCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~DropletSprayDNSThermoCloud() = default;

    // Member Functions
 
        // Access
 
            //- Return a reference to the cloud copy
            inline const thermoCloudType& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;
 
            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();
 
            //- Return the reference value fraction A for mixing
            inline scalar refPropertiesFactor() const;

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

            //- Retruns state of distribution model
            inline const DropletSprayDNSThermo::distributeSourceModel& 
            distributionModel() const;

        // Sources

            //- Return const access to mass source fields
            inline const volScalarField::Internal&
                rhoTrans() const;

            //- Return reference to mass source fields
            inline volScalarField::Internal&
                rhoTrans();

            //- Return const access to enthalpy fields
            inline const volScalarField::Internal&
                hsTrans() const;

            //- Return reference to enthalpy fields
            inline volScalarField::Internal&
                hsTrans();

            //- Return coefficient for carrier phase enthalpy equation
            inline volScalarField::Internal& 
                hsCoeff();

            //- Return const coefficient for carrier phase enthalpy equation
            inline const volScalarField::Internal&
                hsCoeff() const;

            //- Return tmp total mass source for carrier phase
            //- fully explicit
            inline tmp<volScalarField::Internal> Srho() const;

            //- Return total mass source term [kg/m3/s]
            inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;
            
            //- Return total enthalpy term [kg/m3/s]
            inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;

            // Access to moving average function
            inline const movingAverage& averageFunction() const;

        // Provide Source Functions for the DNS field with filtering

            // Note: All source functions take the DNS mesh as an input and 
            // a filtering function which can filter from the mesh of the 
            // cloud to the DNS mesh.

            //- Return tmp total mass source for carrier phase
            //- fully explicit
            inline tmp<volScalarField::Internal> SrhoDNS
            (
                const fvMesh& DNSmesh,
                const filterDNSField& filter
            ) const;

            inline tmp<fvScalarMatrix> SrhoDNS
            (
                const fvMesh& DNSmesh,
                const filterDNSField& filter,
                volScalarField& rho
            ) const;
            
            //- Return total mass fraction Y term [kg/m3/s]
            inline tmp<fvScalarMatrix> SYFuDNS
            (
                const fvMesh& DNSmesh,
                const filterDNSField& filter,
                volScalarField& Yi
            ) const;
            
            //- Return total enthalpy term [kg/m3/s]
            inline tmp<fvScalarMatrix> ShDNS
            (
                const fvMesh& DNSmesh,
                const filterDNSField& filter,
                volScalarField& hs
            ) const;

            inline tmp<fvVectorMatrix> SUDNS
            (
                const fvMesh& DNSmesh,
                const filterDNSField& filter,
                volVectorField& U
            ) const;

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
                const DropletSprayDNSThermoCloud<CloudType,mmcCloudType>& cloudOldTime
            );

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Distribute the sources
            void distributeSources();

            //- Pre-evolve
            void preEvolve(const typename parcelType::trackingData& td);

            //- Evolve the cloud
            void evolve();

        // Integration Scheme Selection

            //- Return the enum type of the integration scheme
            inline const DropletSprayDNSThermo::integrationScheme& 
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
                DropletSprayDNSThermoCloud<CloudType,mmcCloudType>,
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

#include "DropletSprayDNSThermoCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DropletSprayDNSThermoCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
