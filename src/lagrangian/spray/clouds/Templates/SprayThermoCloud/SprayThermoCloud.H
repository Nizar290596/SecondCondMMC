/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SprayThermoCloud

Description
    Templated base class for thermodynamic cloud

    - Adds to kinematic cloud
      - Heat transfer

SourceFiles
    SprayThermoCloudI.H
    SprayThermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef SprayThermoCloud_H
#define SprayThermoCloud_H

#include "thermoCloud.H"
#include "SLGThermo.H"

#include "SprayKinematicCloud.H" 

#include "basicAerosolReactingPopeCloud.H"
#include "interpolationLookUpTable.H"
#include "StochasticLib.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class EvaporationModel;

template<class CloudType>
class CpModel;

template<class CloudType>
class EnvelopeModel;

template<class CloudType>
class FPRadiationModel;

/*---------------------------------------------------------------------------*\
                        Class SprayThermoCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SprayThermoCloud
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
        typedef SprayThermoCloud<CloudType> thermoCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<SprayThermoCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        SprayThermoCloud(const SprayThermoCloud&);

        //- Disallow default bitwise assignment
        void operator=(const SprayThermoCloud&);


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the carrier gas fields

            //- SLG thermodynamics package
            const SLGThermo& thermo_;

            //- Temperature [K]
            const volScalarField& T_;

            //- Pressure [Pa]
            const volScalarField& p_;


        // References to the cloud sub-models

            //-  Evaporation model
            autoPtr<EvaporationModel<SprayThermoCloud<CloudType>>>
                evaporationModel_;

            //- Envelope flame model
            autoPtr<EnvelopeModel<SprayThermoCloud<CloudType>>>
                envelopeModel_;

            //- Cp Model
            autoPtr<CpModel<SprayThermoCloud<CloudType> > >
                cpModel_;

            //- Fuel Particle Radiation Model
            autoPtr<FPRadiationModel<SprayThermoCloud<CloudType> > >
                fpRadiationModel_;
                

        // Reference to the particle integration schemes

            //- Temperature integration
            autoPtr<integrationScheme> TIntegrator_;

            //- Mass integration
            autoPtr<integrationScheme> MIntegrator_;

        // MMC properties data to import

            //-mixture fraction f
            const volScalarField& f_;


        // -LES gradient square of mixture fraction
            const volScalarField& magSqrGradf_;
    
            //- Random mixture fraction 
             StochasticLib1 rand_f_Gen_;
          

            // Reference to the pope cloud
             basicAerosolReactingPopeCloud& pope_;

        //- thermophysical coupling species id table
            List<label> specieId_;

            //- number of liquid fuel Species    
	        label nFuelSpecies_;
            
            //- link between gas and fuel species id
	        scalarField linkFG_;

            //- number of liquid fuel Species    
            label nProductSpecies_;

            //- link between gas and fuel species id
            scalarField linkPG_;

            //- link between gas and fuel species id
            scalarField YProducts_;

        // Sources

            //- Sensible enthalpy transfer [J/kg]
            autoPtr<volVectorField::Internal> MTrans_;

            //- Coefficient for carrier phase hs equation [W/K]
            autoPtr<volScalarField::Internal> MTransRate_;


        //- For treatment of empty supercell event
        List<int> emptySuperCell_;    

        //- For gathering the parcelIds of fp-sp couplings that need rescale timestep
        List<int> parcelId_;

    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(SprayThermoCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        SprayThermoCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& muEff,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            const volScalarField& f,
            basicAerosolReactingPopeCloud& pope,
            const volScalarField& magSqrGradf,
            bool readFields = true
        );

        //- Copy constructor with new name
        SprayThermoCloud(SprayThermoCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        SprayThermoCloud
        (
            const fvMesh& mesh,
            const word& name,
            const SprayThermoCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SprayThermoCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SprayThermoCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~SprayThermoCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const SprayThermoCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return const access to thermo package
            inline const SLGThermo& thermo() const;

            //- Return const access to the carrier temperature field
            inline const volScalarField& T() const;

            //- Return const access to the carrier prressure field
            inline const volScalarField& p() const;

            //- -Return particle mass flow rate
            inline const tmp<volScalarField> mFlow() const;


            // Sub-models

                //- Return reference to heat transfer model
                inline const EvaporationModel<SprayThermoCloud<CloudType>>&
                    evaporation() const;

                //- Return reference to envelope flame model
                inline EnvelopeModel<SprayThermoCloud<CloudType>>&
                    envelope();

                //- Return reference to Cp model
                inline const CpModel<SprayThermoCloud<CloudType> >&
                    cpModel() const;

                //- Return reference to Cp model
                inline const FPRadiationModel<SprayThermoCloud<CloudType> >&
                    fpRadiationModel() const;    

            // Integration schemes

                //-Return reference to velocity integration
                inline const integrationScheme& TIntegrator() const;

                //-Return constant access to mass 
                inline const integrationScheme& MIntegrator() const;


            // -Return const access to the square root of gradient f
            inline const volScalarField& f() const;


            // -LES gradient mixture fraction function
            inline const volScalarField& magSqrGradf() const ;

            // Random mixture fraction
            inline StochasticLib1& rand_f_Gen();      
                  
            //- Return const access to the popeCloud
            inline  basicAerosolReactingPopeCloud& pope() const; 

           
            //- access to species table list id
            inline List<label> specieId() const;

            inline label nFuelSpecies() const;

            inline scalarField linkFG() const;

            inline label nProductSpecies() const;

            inline scalarField linkPG() const;

            inline scalarField YProducts() const;	
            // Sources

                // Enthalpy

                    //- Sensible enthalpy transfer [J/kg]
                    inline volVectorField::Internal& MTrans();

                    //- Sensible enthalpy transfer [J/kg]
                    inline const volVectorField::Internal&
                        MTrans() const;

                    //- Return coefficient for carrier phase hs equation
                    inline volScalarField::Internal& MTransRate();

                    //- Return const coefficient for carrier phase hs equation
                    inline const volScalarField::Internal&
                        MTransRate() const;

                    //- Return tmp momentum source term
                    inline tmp<fvVectorMatrix> SuMass(volVectorField& U) const;

                    //- Return sensible enthalpy source term [J/kg/m3/s]
                    inline tmp<fvScalarMatrix> SM(volScalarField& rho) const;

                    //- Return tmp total mass source for carrier phase  
                    //  - fully explicit   
                    inline tmp<volScalarField::Internal> SM() const;

            //- Return empty supercell event    
            inline List<int> emptySuperCell() const; 

           //- Return constant access of the list of parcelIds
           inline List<int> parcelId() const;

        // Check

            //- Maximum temperature
            inline scalar Tmax() const;

            //- Minimum temperature
            inline scalar Tmin() const;

        // Edit

            //- editing empty sueprcell event
            inline List<int>& emptySuperCell();  

            //- Return editor access to the list of parcelIds 
            inline List<int>& parcelId();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Increase count
            void parcelIdAddOne(const label parcelId); 

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
            void relaxSources(const SprayThermoCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Evolve the cloud
            void evolve();

	    //- set fuel mixture
	    void setFuelMix();

            //- set product mixture
            void setProductMix();

        // Statistics
        
                //- Set Eulerian statistics capturing
                virtual void setEulerianStatistics();

                //- Update Eulerian statistics
                virtual void updateEulerianStatistics();

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

            virtual void writeFields() const;
            
            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayThermoCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SprayThermoCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
