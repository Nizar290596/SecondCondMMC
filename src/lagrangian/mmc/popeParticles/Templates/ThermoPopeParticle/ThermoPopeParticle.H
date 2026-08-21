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
    Foam::ThermoPopeParticle

Description


SourceFiles
    ThermoPopeParticleI.H
    ThermoPopeParticle.C
    ThermoPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoPopeParticle_H
#define ThermoPopeParticle_H

#include "particle.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes 
template<class CloudType, class ReactionThermo>
class ThermoPopeCloud;

template<class ParticleType>
class ThermoPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const ThermoPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class ThermoPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ThermoPopeParticle
:
    public ParticleType
{
    // Friend declarations
        template<class CloudType, class ReactionThermo> 
        friend class ThermoPopeCloud;  

public:
    
    typedef ThermoPopeParticle<ParticleType> thermoParticleIOType;
    
    class trackingData
    :
        public ParticleType::trackingData
    {
        private:

        // Private data

            // Interpolators for continuous phase fields

                //- Pressure field interpolator
                autoPtr<interpolation<scalar> > pInterp_;

        public:

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                TrackCloudType& cloud
            );


        // Member functions

            //- Return const access to the interpolator for continuous
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;
    };
    
    
    // Public member data
        //- Ratio of specific heats[]
        static constexpr scalar gamma = 1.4;

protected:

    // Protected data

       //Static data
            static HashTable<label, word> indexInXiC_;
            
            static List<word> XiCNames_;
            
            static List<word> componentNames_;

        // Particle properties

            //- Temperature [K]
            scalar T_;
            
            //- Mass fractions of mixture []
            scalarField Y_;
            
            // - Coupling Variables
            scalarField XiC_;
            
            //- Total enthalpy [J/kg]
            scalar hA_;
  
            //- Equivalent enthalpy [J/kg]
            scalar hEqv_;
            
            //- Species source term from the dispersed phase []
            scalarField dpYsource_;

            //- Mixture fraction source term from the dispersed phase []
            scalarField dpXiCsource_;

            //- Enthalpy source term from the dispersed phase [J/kg]
            scalar dp_hAsource_;

        // Cell-based quantities
        
            //- Pressure [Pa]
            scalar pc_;

            //- number of active species []
            label NumActSp_;

public:

    // Static data members

        //- Runtime type information
        TypeName("ThermoPopeParticle");

        //- String representation of properties
        AddToPropertyList
        (
            ParticleType,
            " T"
          + " Y"
          + " XiC"
          + " hA"
          + " hEqv"
          + " dpYsource"
          + " dpXiCsource"
          + " dp_hAsource"
          + " pc"
          + " NumActSp"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ThermoPopeParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        ThermoPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        ThermoPopeParticle(const ThermoPopeParticle& p);

        //- Construct as a copy
        ThermoPopeParticle(const ThermoPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ThermoPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ThermoPopeParticle(*this, mesh));
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

            autoPtr<ThermoPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<ThermoPopeParticle<ParticleType> >
                (
                    new ThermoPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

	    //- Second-conditioning flag stub.
            //  Real value lives on MixingPopeParticle; this base-layer stub
            //  lets generic templates (e.g. ThermoPhysicalCouplingModel)
            //  reference it on particle types that do not carry a subset
            //  flag. Always returns 0 (not in the subset).
            inline label secondCondFlag() const { return 0; }

	    //- Modified progress variable (phi-degree) stub.
            //  Real value lives on MixingPopeParticle; this base-layer stub
            //  lets generic templates (e.g. ThermoPhysicalCouplingModel)
            //  reference it on particle types that do not carry it.
            inline scalar phiModified() const { return 0.0; }

            //- Return const access to total enthalpy
            inline const scalar& hA() const;
            
            //- Return const access to target enthalpy
            inline const scalar& hEqv() const;

            //- Return const access to mass fractions of mixture
            inline const scalarField& Y() const;

            //- Return const access to species source term from dispersed phase
            inline const scalarField& dpYsource() const;

            //- Return const access to enthalpy source term from dispersed phase
            inline const scalar& dp_hAsource() const;

            //- Return const access to coupling variables
            inline const scalarField& XiC() const;
            
            //- Return const access to specific coupling variable
            inline const scalar& XiC(const word& name) const;
            
            //- Return const access to mixture fraction source term from 
            //- dispersed phase
            inline const scalarField& dpXiCsource() const;
            
            //- Return const access to specific coupling variable source term 
            //- from dispersed phase
            inline const scalar& dpXiCsource(const word& name) const;

            //- Return the owner cell pressure
            inline scalar pc() const;   

            //- Return const access to number of species
            inline const label& NumActSp() const;



        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return const access to total enthalpy
            inline scalar& hA();
            
            //- Return const access to target enthalpy
            inline scalar& hEqv();

            //- Return access to mass fractions of mixture
            inline scalarField& Y();

            //- Return access to species source term from dispersed phase
            inline scalarField& dpYsource();

            //- Return access to enthalpy source term from dispersed phase
            inline scalar& dp_hAsource();
            
            //- Return access to coupling variables
            inline scalarField& XiC();
            
            //- Return const access to specific coupling variable
            inline scalar& XiC(const word& name);
            
            //- Return access to mixture fraction source term from dispersed phase
            inline scalarField& dpXiCsource();
            
            //- Return const access to specific coupling variable source term 
            //- from dispersed phase
            inline scalar& dpXiCsource(const word& name);

            //- Return reference to the owner cell pressure
            inline scalar& pc();

            //- Return access to number of species
            inline label& NumActSp();
            
            template<class CloudType>
            static void setStaticProperties(CloudType& c); 


        // Main calculation loop

            //- Update particle properties from the finite volume field
            template<class TrackCloudType>
            void setCellValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt,
                const label cellI
            );
            
            //- Update particle properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt,
                const label cellI
            );

            //- Update particle properties from the dispersed phase
            inline void updateDispersedSources();
            
            //- Mix operator for this particleType
            static void mixProperties
            (
                ThermoPopeParticle<ParticleType>& p, 
                ThermoPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            static void mixProperties
            (
                ThermoPopeParticle<ParticleType>& p, 
                ThermoPopeParticle<ParticleType>& q, 
                scalar mixExtent,
                scalarList ScaledExtent
            );
        
        // For particle sampling 

            //- Initialize settings for statistical data sampling
            void initStatisticalSampling();

            //- Return the data in a container
            DynamicList<scalar> getStatisticalData(const wordList& vars) const;

            //- Return variables names stored in the container
            DynamicList<word> getStatisticalDataNames(const wordList& vars) const;

        // I-O

            //- Read particle data for all particles in a cloud from files
            //  Use given composition model
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read particle data for all particles in a cloud from files
            //  Without composition model
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write particle data for all particles in a cloud from files
            //  Use given composition model
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write particle data for all particles in a cloud from files
            //  without composition model
            template<class CloudType>
            static void writeFields(const CloudType& c);

    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const ThermoPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoPopeParticleI.H"
#include "ThermoPopeParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

