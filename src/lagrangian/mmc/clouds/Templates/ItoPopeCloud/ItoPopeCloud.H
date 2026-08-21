/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ItoPopeCloud

Description
    Contains the ItoPopeParticles, with the stochastic movement.
    Has the move() function which calls calc() on the particles to calculate
    the chemical source term. 

    Can swith on massConsistencyCheck which writs out the mass balance for 
    each time step into the file 

SourceFiles
    ItoPopeCloudI.H
    ItoPopeCloud.C

Author
Edison GE
y.ge1222@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef ItoPopeCloud_H
#define ItoPopeCloud_H

#include "particle.H"
#include "Cloud.H"
#include "itoPopeCloud.H"
#include "IOdictionary.H"
#include "StochasticLib.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "scalarIOField.H"
#include "popeCloudSolution.H"
#include "eulerianStatistics.H"

#include "turbulentFluidThermoModel.H"
#include "particleNumberController.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class InflowBoundaryModel;

 
/*---------------------------------------------------------------------------*\
                        Class ItoPopeCloud Declaration
\*---------------------------------------------------------------------------*/
template<class CloudType>
class ItoPopeCloud
:
    public CloudType,
    public itoPopeCloud
{
    public:
        
    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of particle the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef ItoPopeCloud<CloudType> itoPopeCloudType;

    private:
    
        //- Cloud copy pointer
        autoPtr<ItoPopeCloud<CloudType> > cloudCopyPtr_;

    protected:
        
        //- References to the mesh
        const fvMesh& mesh_;
        
        // Dictionaries for settings 
        
            //- Dictionary of cloud properties
            IOdictionary cloudProperties_;
            
            //- Solution properties
            popeCloudSolution solution_;
                    
            //- Sub-models dictionary
            const dictionary subModelProperties_;

        // Flags

            //- Switch to write out the mass consistency check 
            Switch writeMassConsistency_;
        
        // Particle number management
        
            //- Particle management class
            particleNumberController pManager_;
            
        // Other variables 
        
            //- Reference to time and scalar data for time ?? at previous timestep ??
            const Time& runTime_;
            scalar time0_;

            //- Reference to the turbulence model (LES)
            const compressible::LESModel* turbModelPtr_;

            //- Random numbers
            StochasticLib1 rndGen_;

            //- In-cloud cache of deltaT, lookup in submodels and parcel is expensive
            scalar cachedDeltaT_;

        //- References to the physical fields
        
            //- Eulerian velocity
            const volVectorField& U_;
            
            //- Eulerian diffusivity
            const volScalarField& DEff_;
            
            //- Eulerian density
            const volScalarField& rho_;
            
            //- Eulerian grad(rho_*DEff)
            //  Required for the mixing models and the drift term
            const volVectorField& gradRho_;

        //- Cell length scale
        scalarField cellLengthScale_;

        //- Standard mass represented by Pope particle with weight wt=1.0.
        //  This is just for normalisation and any mass can be used.
        static constexpr scalar deltaM_ = 1.0;
        
        //- Eulerian statistics
        const dictionary eulerianStatsDict_;
        eulerianStatistics eulerianStats_;
        
        // Reference to the cloud sub-models
        //- Inflow boundary model
        autoPtr<InflowBoundaryModel<ItoPopeCloud<CloudType> > > inflowBoundaryModel_;

        //- Finite volume to particle mass difference
        scalar fvMinusPpMass_;

        //- File stream to write out the mass consistency check
        autoPtr<OFstream> massConsistencyCheck_;

    // Private member functions

        //- Lookup field or construct
        template<class Type>
        GeometricField<Type, fvPatchField, volMesh>&
        lookupOrConstruct(const word fieldName, dimensionSet units=dimless);

        //- Construct cloud dictionary name from cloudName
        //  the cloudName is by default: popeParticlesSet_1, popeParticlesSet_2,...
        //  which results then in the dicionaries:
        //  popeParticleSet_1 : cloudProperties
        //  popeParticleSet_2 : cloudProperties_2
        //  popeParticleSet_3 : cloudProperties_3
        //  the first cloudProperties is kept for backwards compatibility
        word createCloudPropertiesName(const word&) const;

    public:

    // Constructors

        //- Construct with given fields
        ItoPopeCloud
        ( 
            const word& cloudName,
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& DEff,
            const volScalarField& rho,
            const volVectorField& gradRho,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

        //- Construct and lookup fields in object registry
        ItoPopeCloud
        ( 
            const word& cloudName,
            const fvMesh& mesh,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

        //- Disallow default bitwise copy construct
        ItoPopeCloud(const ItoPopeCloud&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ItoPopeCloud&) = delete;
     
     
        // Destructor
        virtual ~ItoPopeCloud();


        // Member Functions

            // Initialisation

                //- Set cloud sub-models
                void setModels();

            // Cloud evolution functions

                //- Solve the cloud - calls all evolution functions
                template<class TrackCloudType> 
                void solve
                (
                    TrackCloudType& cloud,
                    typename particleType::trackingData& td
                );

                //- Reset state of cloud
                void cloudReset(ItoPopeCloud<CloudType>& c);

                //- Initial release of particles
                void initReleaseParticles();

                //- Move the particles
                template<class TrackCloudType>
                void move
                (
                    TrackCloudType& cloud,
                    typename particleType::trackingData &td
                );
            
            //- Calculate fv and particle mass
            void checkMassConsistency
            (
                scalar& fvM,
                scalar& ppM
            );
 
            inline scalar& fvMinusPpMass();
            
            inline void storeDeltaT();
            inline scalar cachedDeltaT() const;

            // Accessors
            
            //sub-models
                //- Return constant reference to inflow boundary
                inline const InflowBoundaryModel<ItoPopeCloud<CloudType> >& inflowBoundary() const;

                //- Return non-constant reference to inflow boundary
                inline InflowBoundaryModel<ItoPopeCloud<CloudType> >& inflowBoundary();
            
            
            inline const fvMesh& mesh() const;
            inline scalar deltaM() const;
            inline const Time& runTime() const;
            inline StochasticLib1& rndGen();
            
            //- Return flag to write out the massConsistencyCheck
            const Switch& writeMassConsistency() const;


            //- References to the carrier gas fields

            inline const volVectorField& U() const;
            inline const volScalarField& DEff() const;
            inline const volScalarField& rho() const;
            inline const volVectorField& gradRho() const;
            inline const volScalarField& delta() const;
            inline const volScalarField& muEff() const;
            inline const volScalarField& kfv() const;//!! turbulent kinetic energy

            //- Return dictionary name
            inline const IOdictionary& cloudProperties() const;
            inline const dictionary& subModelProperties() const;
            inline popeCloudSolution& solution();
            inline const popeCloudSolution& solution() const;

            //- Return cell length
            inline const scalarField& cellLengthScale() const;

            //- Return statistics
            inline eulerianStatistics& eulerianStats();
            inline const eulerianStatistics& eulerianStats() const;
            
            //- Return statistics dictionary
            inline const dictionary& eulerianStatsDict() const;

            //- Clear the cloud
            inline void clear();
            
            //-  Add new parcel
            void addNewParticle
            (
                const barycentric& coordinates,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const scalar wt,
                const scalar mass,
                const label patchI=-1,
                const label patchFace=-1,
                const bool iniRls=true
            );
            
            //- Set particle thermo properties
            virtual void setParticleProperties
            (
                particleType& particle,
                const scalar& mass,
                const scalar& wt,
                const scalar& patchI,
                const scalar& patchFace,
                const bool& iniRls
            );
            
            //- Set Eulerian statistics capturing
            virtual void setEulerianStatistics();

            //- Update Eulerian statistics
            virtual void updateEulerianStatistics();
            
            virtual void writeFields() const;
            
        // Return particle manager
        const particleNumberController& pManager() const {return pManager_;}
};

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ItoPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ItoPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
