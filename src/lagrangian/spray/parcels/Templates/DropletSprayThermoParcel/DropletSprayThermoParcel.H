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
    Foam::DropletSprayThermoParcel

Description
    Thermodynamic parcel class representing single component droplets in a 
    spray. Includes coupling to the stochastic particles and considers 
    evaporation.

SourceFiles
    DropletSprayThermoParcelI.H
    DropletSprayThermoParcel.C
    DropletSprayThermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef DropletSprayThermoParcel_H
#define DropletSprayThermoParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"
#include "nameVariableTable.H"

// We need to know about the enum type
#include "DropletSprayThermoCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class DropletSprayThermoParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const DropletSprayThermoParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                    Class DropletSprayThermoParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DropletSprayThermoParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


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

            //- Particle specific heat capacity [J/(kg.K)]
            demandDrivenEntry<scalar> Cp0_;

            //- Particle emissivity [] (radiation)
            demandDrivenEntry<scalar> epsilon0_;

            //- Particle scattering factor [] (radiation)
            demandDrivenEntry<scalar> f0_;

    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member functions

            // Access

                //- Return const access to the particle initial temperature [K]
                inline scalar T0() const;

                //- Return const access to minimum temperature [K]
                inline scalar TMin() const;

                //- Return const access to maximum temperature [K]
                inline scalar TMax() const;

                //- Set the maximum temperature [K]
                inline void setTMax(const scalar TMax);

                //- Return const access to the particle specific heat capacity
                //  [J/(kg.K)]
                inline scalar Cp0() const;

                //- Return const access to the particle emissivity []
                //  Active for radiation only
                inline scalar epsilon0() const;

                //- Return const access to the particle scattering factor []
                //  Active for radiation only
                inline scalar f0() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            //- Local copy of carrier specific heat field
            //  Cp not stored on carrier thermo, but returned as tmp<...>
            const volScalarField Cp_;

            //- Local copy of carrier thermal conductivity field
            //  kappa not stored on carrier thermo, but returned as tmp<...>
            const volScalarField kappa_;


            // Interpolators for continuous phase fields

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;

                //- Reference variable list of field interpolators
                HashTable
                <
                    autoPtr<interpolation<scalar>>
                >
                XiRInterp_;

                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpInterp_;

                //- Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaInterp_;

                //- Radiation field interpolator
                autoPtr<interpolation<scalar>> GInterp_;


            // Cached continuous phase properties

                //- Temperature [K]
                scalar Tc_;

                //- Specific heat capacity [J/(kg.K)]
                scalar Cpc_;

    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions

            //- Return access to the locally stored carrier Cp field
            inline const volScalarField& Cp() const;

            //- Return access to the locally stored carrier kappa field
            inline const volScalarField& kappa() const;

            //- Return const access to the interpolator for continuous
            //  phase temperature field
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for named continuous
            //  reference variables
            inline const interpolation<scalar>& XiRInterp
            (
                const word name
            ) const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase thermal conductivity field
            inline const interpolation<scalar>& kappaInterp() const;

            //- Return const access to the interpolator for continuous
            //  radiation field
            inline const interpolation<scalar>& GInterp() const;

            //- Return the continuous phase temperature
            inline scalar Tc() const;

            //- Access the continuous phase temperature
            inline scalar& Tc();

            //- Return the continuous phase specific heat capacity
            inline scalar Cpc() const;

            //- Access the continuous phase specific heat capacity
            inline scalar& Cpc();
    };


private:

    // Private Member Functions

        //- Transform a species mass fraction vector Y to the mole fraction X
        List<scalar> YtoX(const List<scalar>& Y, const List<scalar>& MW) const;

        //- Transform a species mass fraction vector Y to the mole fraction X
        List<scalar> XtoY(const List<scalar>& X, const List<scalar>& MW) const;

        //- Integrate the evaporation in time with an Euler explicit method
        template<class TrackCloudType>
        void timeIntegrationEulerExplicit
        (
            TrackCloudType& cloud,
            trackingData& td,
            scalar& muRef,              // Viscosity at reference conditions
            scalar& rhoRef,             // Gas density at reference conditions
            scalar& mDNew,              // Mass of the droplet after evaporation
            scalar& TDNew,              // Droplet temperature after evaporation
            const scalar pInf,          // Infinity pressure
            const scalar Ts,            // Surface temperature
            const scalar TInf,          // Infinity temperature
            const List<scalar>& YInf,   // Infinity species composition
            const scalar dt             // Time step
        ) const;


        //- Integrate the evaporation in time with Crank-Nicolson method
        template<class TrackCloudType>
        void timeIntegrationCrankNicolson
        (
            TrackCloudType& cloud,
            trackingData& td,
            scalar& muRef,              // Viscosity at reference conditions
            scalar& rhoRef,             // Gas density at reference conditions
            scalar& mDNew,              // Mass of the droplet after evaporation
            scalar& TDNew,              // Droplet temperature after evaporation
            const scalar pInf,          // Infinity pressure
            const scalar Ts,            // Surface temperature
            const scalar TInf,          // Infinity temperature
            const List<scalar>& YInf,   // Infinity species composition
            const scalar dt             // Time step
        ) const;

        // Solves a 2x2 scalar matrix system Ax = b
        Pair<scalar> matrixSolver2By2
        (
            const scalar a11,
            const scalar a12,
            const scalar a21,
            const scalar a22,
            const scalar b1,
            const scalar b2
        ) const;

protected:

    // Protected data

        // For particle statistics

            //- Allow name variable access -- used for particle sampling
            nameVariableTable table_{"DropletSprayParcel"};

        // Parcel properties

            //- Temperature [K]
            scalar T_;

            //- Specific heat capacity [J/(kg.K)]
            scalar Cp_;

            //- Mass of the droplet
            scalar mass0_;

            //- Species index in the composition
            //  Default initialize to -1
            label jFuel_{-1};

            //- Evaporation mass flux
            //  stored for particle statistics
            scalar mDot_{0};

public:

    // Static data members

        //- Runtime type information
        TypeName("DropletSprayThermoParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " T"
          + " Cp"
          + " mass0"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline DropletSprayThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline DropletSprayThermoParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline DropletSprayThermoParcel
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
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        DropletSprayThermoParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        DropletSprayThermoParcel(const DropletSprayThermoParcel& p);

        //- Construct as a copy
        DropletSprayThermoParcel(const DropletSprayThermoParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new DropletSprayThermoParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new DropletSprayThermoParcel(*this, mesh));
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

            autoPtr<DropletSprayThermoParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<DropletSprayThermoParcel<ParcelType>>
                (
                    new DropletSprayThermoParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

            //- Return const access to specific heat capacity
            inline scalar Cp() const;

            //- Return the parcel sensible enthalpy
            inline scalar hs() const;

            //- Return access to species index in Y
            inline label jFuel() const;

            //- Return access to the droplet mass stored on the particle
            //  Note: this->mass() returns the mass calculated based on the volume
            inline scalar mass0() const;


        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return access to specific heat capacity
            inline scalar& Cp();

            //- Return access to the droplet mass stored on the particle
            //  Note: this->mass() returns the mass calculated based on the volume
            inline scalar& mass0();

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Calculate reference properties based on the 1/3 rule
            template<class TrackCloudType>
            void calcReferenceValues
            (
                TrackCloudType& cloud,
                const scalar pInf,          // Infinity pressure
                const scalar Ts,            // Surface temperature
                const scalar TInf,          // Infinity temperature
                const List<scalar>& Ys,     // Surface species composition
                const List<scalar>& YInf,   // Infinity species composition
                const label jFuel,          // Index of the fuel mass fraction in Y
                scalar& TRef,               // Reference temperature
                scalar& rhoRef,             // Density at reference conditions
                scalar& muRef,              // Viscosity at reference conditions
                scalar& kappaRef,           // heat conductivity at reference condition
                scalar& CpRef,              // heat capacity at reference conditions
                scalar& CpFRef              // heat capacity at reference conditions (pure liquid mass fraction)
            ) const;

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Calculate evaporation rate [kg/s]
            template<class TrackCloudType>
            void mDotAndQDot
            (
                TrackCloudType& cloud,
                trackingData& td,
                scalar& mDotF,              // Evaporated mass of fuel [kg/s]
                scalar& QDot,               // Heat flux [J/s]
                scalar& muRef,              // Viscosity at reference conditions
                scalar& rhoRef,             // Gas density at reference conditions
                const scalar dD,            // Droplet diameter
                const scalar pInf,          // Infinity pressure
                const scalar Ts,            // Surface temperature
                const scalar TInf,          // Infinity temperature
                const List<scalar>& YInf,   // Infinity species composition
                const label jFuel           // Index of the fuel mass fraction in Y
            ) const;

        // Check

            //- Resets the mDot on the particle
            inline void reset();

        // For Particle Sampling

            //- Access to the name-variable table 
            nameVariableTable& nameVariableLookUpTable() {return table_;};

            const nameVariableTable& nameVariableLookUpTable() const {return table_;};

            //- Initialize settings for statistical data sampling
            void initStatisticalSampling();

            //- Return the data in a container
            //  Optionally provide a list of variables to sample
            DynamicList<scalar> getStatisticalData(const wordList& vars) const;

            //- Return variables names stored in the container
            //  Optionally provide a list of variables to sample
            DynamicList<word> getStatisticalDataNames(const wordList& vars) const;

        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
            ) const;

            //- Read particle fields as objects from the obr registry
            template<class CloudType>
            static void readObjects(CloudType& c, const objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            template<class CloudType>
            static void writeObjects(const CloudType& c, objectRegistry& obr);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const DropletSprayThermoParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DropletSprayThermoParcelI.H"
#include "DropletSprayThermoParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DropletSprayThermoParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
