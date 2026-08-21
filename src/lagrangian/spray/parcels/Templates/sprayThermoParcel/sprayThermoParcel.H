/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::sprayThermoParcel

Description
    Thermodynamic parcel for sprays to add 


SourceFiles
    sprayThermoParcelI.H
    sprayThermoParcel.C
    sprayThermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef sprayThermoParcel_H
#define sprayThermoParcel_H

#include "fieldTypes.H"
#include "particle.H"
#include "interpolationLookUpTable.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"
#include "interpolation.H"
#include "IOstream.H"
#include "basicAerosolReactingPopeCloud.H"
#include "nameVariableTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class sprayThermoParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const sprayThermoParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class sprayThermoParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class sprayThermoParcel
:
    public ParcelType
{
    // Private data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold thermo particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private data

            //- Particle initial temperature [K]
            demandDrivenEntry<scalar> T0_;

            //- Minimum temperature [K]
            demandDrivenEntry<scalar> TMin_;

            //- Maximum temperature [K]
            demandDrivenEntry<scalar> TMax_;

            //-Schmidt number
            demandDrivenEntry<scalar> Sc_;

            //-Lewis number
            demandDrivenEntry<scalar> Le_;

            //- latent heat
            demandDrivenEntry<scalar> LH_;

            //- formation heat
            demandDrivenEntry<scalar> FH_;
            
            //-Molecular weight of surface state specei
            demandDrivenEntry<scalar> molWtSstate0_;
        
    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties
            (
                const dictionary& parentDict,
                const bool readFields = true
            );

            //- Construct from components
            constantProperties
            (
                const label parcelTypeId,
                const scalar rhoMin,
                const scalar minParticleMass,
                const scalar T0,
                const scalar TMin,
                const scalar TMax,
                const scalar Sc,
                const scalar Le,
                const scalar latentHeat,
                const scalar formationHeat,
                const scalar molWtSstate0
	    );


        // Member functions

            // Access

                //- Return const access to the particle initial temperature [K]
                inline scalar T0() const;

                //- Return const access to minimum temperature [K]
                inline scalar TMin() const;

                //- Return const access to maximum temperature [K]
                inline scalar TMax() const;
          
                //- Return const access to the default carrier schmidt number []
                inline scalar Sc() const;

                //- Return const access to the default carrier Lewis number []
                inline scalar Le() const;

                //- Return const access to the default latent heat
                inline scalar LH() const;

                //- Return const access to the default latent heat
                inline scalar FH() const;
                
                //- Return const access to the default molecular weight of surface specei
                inline scalar molWtSstate0() const;

    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;

                //- mixtrue fraction field interpolator
                autoPtr<interpolation<scalar>> fInterp_;

                //- pressure field interpolator
                autoPtr<interpolation<scalar>> PInterp_;

                //- Gradient of mixtrue fraction field interpolator
                autoPtr<interpolation<scalar>> sqrGradInterp_;


    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions

            //- phase temperature field
            inline const interpolation<scalar>& TInterp() const;
            
            // Return const access to the interpolator of mixture fraction   
            inline const interpolation<scalar>& fInterp() const; 

            // Return access to Ambient pressure at the location of the fuel particle 
            inline const interpolation<scalar>& PInterp() const; 

            // Return const access to the interpolator of mixture fraction   
            inline const interpolation<scalar>& sqrGradInterp() const;
    };


protected:

    // Protected data

            //- Allow name variable access -- used for particle sampling
            nameVariableTable table_;

            //- Temperature [K]
            scalar TD_;

            //- Specific Heat of the Dispersed Phase[KJ/Kg K]
            scalar Cp_;

            //- Correct Temperature of Gas
            scalar TGas_;

            //- Temp of fuel particle locations
            scalar TFilm_;
 
            scalar TSurf_;

        // Cell-based quantities

            //- pressue of fuel particle surface (pa)
            scalar pc_;

            //- Accurate mixture fraction at the location of fuel particle
            scalar fGas_;
            
            //- Accurate mixture fraction at surface 
            scalar fSurf_;

            //- Accurate mixture fraction at film state
            scalar fFilm_;
        
            //- Number of species
            label nThermoSpecies_;

            //- Number of fuel species
            label nFuelSpecies_;

	        //- Fuel mass fraction fields        
            scalarField YFuel_;

            //- Fuel vapor mass fraction fields        
            scalarField YFuelVap_;

            //- Gas mass fraction fields        
            scalarField YGas_;

            //- Filmstate specie
            scalarField YFilm_;

            //-Surface specieI
            scalarField YSurf_;             

            //-standardized Gas state enthalpy
            scalar haGas_;

            //-standardized dispersion state enthalpy
            scalar hD_;
        
            //-standardized surface state enthalpy
            scalar haSurf_;

            //-Standardized film enthalpy
            scalar haFilm_;
            
            //-Standardized enthalpy of the gaseous fuel composition
            scalar haFuelVap_;
            
            //- Species based viscoisty
            scalar muFilm_;

            //- Species based density
            scalar rhoFilm_;
            
            //-Internal heating
            scalar qR_;
           
            //-Internal heating
            scalar qD_;
            
            //-massFlux
            scalar mFlux_;
    
            //-Molecular wt of surface specei
            scalar molWtSstate_;  

            //- film Z
            scalar ZFilm_;

            //- spalding number
            scalar BM_;

            scalar superCelli_;

            //- switch for saturated gas phase
            label noMassTransfer_;

public:

    // Static data members

        //- Runtime type information
        TypeName("sprayThermoParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
              " TD "             
            + " TGas " 
            + " Cp " 
            + " TFilm "
            + " TSurf "
            + " pc "
            + " fGas "
            + " fSurf "
            + " fFilm "
            + " nThermoSpecies"
            + " nFuelSpecies"
            + " YFuel"
	        + " YFuelVap"
            + " YGas "
            + " YFilm "
            + " YSurf "
            + " haGas "
            + " hD "  
            + " haSurf "
            + " haFilm "
            + " haFuelVap " 
            + " muFilm "
            + " rhoFilm "
            + " YM "
            + " qD "
            + " mFlux "
            + "molWtSstate "
            + "ZFilm"
            + "BM"
        );

    // Public type defs
        //- Type of particle the pope cloud was instanciated for 
        typedef typename basicAerosolReactingPopeCloud::particleType scParticleType;

    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline sprayThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline sprayThermoParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline sprayThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        sprayThermoParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        sprayThermoParcel(const sprayThermoParcel& p);

        //- Construct as a copy
        sprayThermoParcel(const sprayThermoParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new sprayThermoParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new sprayThermoParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<sprayThermoParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<sprayThermoParcel<ParcelType>>
                (
                    new sprayThermoParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

            //- Return const access to Specific Heat of Dispersed Phase
            inline scalar Cp() const;

            //- Return the parcel sensible enthalpy
            inline scalar hs() const;

            //- Return const access to Gas temperature
            inline scalar TGas() const;

            //- Return const access to surface Temperature of Dispersed Phase
            inline scalar TSurf() const;

            //- Return const access to film temperature
            inline scalar TFilm() const;

            //- Return const access to carrier specific heat capacity
            inline scalar Cpc() const;
            
            //-Return const access to pressure
            inline scalar pc() const;

            //- Return Accurate mixture fraction at the location of fuel particle
            inline scalar fGas() const ;

            //- Return Accurate mixture fraction at the location of fuel particle
            inline scalar fSurf() const ;

            //- Return const access to filmstate mixtrue fraction
            inline scalar fFilm() const;

            //- Access to thermospecies
            inline label nThermoSpecies() const;

            //- Access to thermospecies
            inline label nFuelSpecies() const;

            //- Return const access to species field  
            inline scalarField YFuel() const;

            //- Return const access to species field  
            inline scalarField YFuelVap() const;

            //- Return const access to species field  
            inline scalarField YGas() const;

            //- Return const access to species mass fraction
            inline  scalarField YFilm() const;    

            //- Return const access to species mass fraction
            inline  scalarField YSurf() const;    

            //- Return const access to standardized enthalpy
            inline scalar haGas() const; 

            //- Return const access to standardized enthalpy for dispersion state
            inline scalar hD() const;
            //- Return const access to standardized enthalpy
            inline scalar haSurf() const; 
      
            //- Return const access to standardized enthalpy of film state 
            inline scalar haFilm() const;
        
            //- Return const access to standardized enthalpy of middle state
            inline scalar haFuelVap() const;

            //- Return const access to film viscoisty
            inline scalar muFilm() const;
           
            //- Return const access to film density
            inline scalar rhoFilm() const; 

            //-Return const access to mass fraction
            inline scalarField YM() const;   

            //-Internal heating
            inline scalar qR() const;

            //-Internal heating
            inline scalar qD() const;    
            
            //-Access to mass flux
            inline scalar mFlux() const;            

            //-Access to mol wt of surface speciei
            inline scalar molWtSstate() const;

             //- Access to mass transfer rate
            inline scalar ZFilm() const;       

             //- Access to mass transfer rate
            inline scalar BM() const;  

        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return  access to Specific Heat of Dispersed Phase
            inline scalar& Cp() ;

            inline scalar& TGas();

            //- Return access to number of thermo species
            inline label& nThermoSpecies();

            //- Return access to number of fuel  species
            inline label& nFuelSpecies();

            inline scalarField& YFuel();

            inline scalar& molWtSstate() ;

            //-Calculate Nusselt Number
            inline scalar NuFilm
            (   
                const scalar ReFilm,
                const scalar PrFilm
            ) const;

            //-Calculate Sherwood Number           
            inline scalar ShFilm
            (   
                const scalar ReFilm,
                const scalar Sc
            ) const;
             
            //-Calculate Prandtl Number 
            inline scalar PrFilm
            (   
                const scalar muFilm,
                const scalar rhoFilm,
                const scalar alpha
            ) const;   

            //- Calculate Spalding number for mass transfer
            inline scalar BM
            (
                const scalar fGas,           // Gas mixture fraction
                const scalar fSurf           // Surface surface mixture fraction
            ) const;

            //- Calculate Spalding number for mass transfer
            inline scalar BH
            (
                const scalar haGas,           // Gas enthalpy
                const scalar haSurf,          // Surface enthalpy
                const scalar hD,              //-Dispersion phase enthalpy 
                const scalar qD,              //-Internal heating  
                const scalar qR              //- heat due to radiation   
            ) const;
            
            inline scalar temperatureUpdateSwitch(scalar F0, scalar F1);

            //- mmc source term 
            inline scalar calcMMCSource
            (
                const scalar dt, 
                const scalar c, 
                const scalar phi_n, 
                const scalar tau
            );  
   
            //- Calculate radiative heat transfer
            inline scalar qR
            (
                const scalar TD,
                const scalar TWall,
                const scalar mFlux
            ) const;

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues
            (
                TrackCloudType& cloud,
                trackingData& td
            );

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            // Calculate thermophysical properties
            template<class TrackCloudType>
            void calcThermophysicalValues
            (
                TrackCloudType& cloud,
                trackingData& td
            );

            // Calculate surface mixture fraction
            template<class TrackCloudType>
            void calcSurfaceProperties
            (
                TrackCloudType& cloud,
                trackingData& td 
            );

	        template<class TrackCloudType>
            void particleH
            (
                TrackCloudType& cloud,
                trackingData& td
            );

            //- turns compsition of n_fuel components into n_thermo components
            template<class TrackCloudType>
	        scalarField YLiqToYGas
	        (
                TrackCloudType& cloud,
	            trackingData& td,
	            scalarField YLiq
            );

            template<class TrackCloudType>
            void calcSourceLES
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt,
                vector& dMTrans,            
                scalar& dMTransRate,                
                const label celli
            );  
    
            //- Calculate new particle temperature and diameter
            template<class TrackCloudType>
            void calcParticleUpdate
            (
                TrackCloudType& cloud,
                trackingData& td,  
                scalar M0, 
                scalar F0, 
                scalar msp, 
                scalar& hGTG, 
                scalar& hGTD, 
                scalar& hDvol,
                scalar& dm_fp,
                scalar& timeFactor,
                const scalar dt  
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // For particle sampling 

            //- Access to the name-variable table 
            nameVariableTable& nameVariableLookUpTable() {return table_;}

            //- Const access to the name-variable table 
            const nameVariableTable& nameVariableLookUpTable() const {return table_;}

            //- Initialize settings for statistical data sampling
            void initStatisticalSampling();

            //- Return the data in a container
            DynamicList<scalar> getStatisticalData(const wordList& vars) const;

            //- Return variables names stored in the container
            DynamicList<word> getStatisticalDataNames(const wordList& vars) const;

        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const sprayThermoParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sprayThermoParcelI.H"
#include "sprayThermoParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sprayThermoParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
