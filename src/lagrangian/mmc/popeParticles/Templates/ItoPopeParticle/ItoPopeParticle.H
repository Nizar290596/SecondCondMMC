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
    Foam::ItoPopeParticle

Description
    Kinematic parcel class with rotational motion (as spherical
    particles only) and one/two-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions

SourceFiles
    ItoPopeParticleI.H
    ItoPopeParticle.C
    ItoPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef ItoPopeParticle_H
#define ItoPopeParticle_H

#include "particle.H"
#include "IOstream.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"
#include "autoPtr.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "Random.H"
#include "contiguous.H"
#include "baseParticleDataContainer.H"
#include "nameVariableTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class ItoPopeParticle;

// Forward declaration of friend functions

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const ItoPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                         Class ItoPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ItoPopeParticle
:
    public ParticleType
{
public:
    
    typedef ItoPopeParticle<ParticleType> itoParticleIOType;
    
    class trackingData
    :
        public ParticleType::trackingData
    {
        //- Counters for the statistics
        label wallCollisions_;
        label wedgeCollisions_ ;
        label symmetryPlaneCollisions_ ;
        label leavingModel_;
        label injectedInModel_;
        label changedProzessor_;
        label noOutflowRebound_;
        label noRandomOutflowRebound_;
        label periodicModel_;

        // Interpolators for continuous phase fields

        autoPtr<interpolation<vector> >  UInterp_;
        autoPtr<interpolation<scalar> >  DInterp_;
        autoPtr<interpolation<scalar> >  rhoInterp_;
        autoPtr<interpolation<vector> >  gradRhoInterp_;

    public:

        //- Constructor
        template<class TrackCloudType>
        trackingData
        (
            const TrackCloudType& cloud
        );

        //- Access

        void resetCounters()
        {
            wallCollisions_=0;
            wedgeCollisions_=0;
            symmetryPlaneCollisions_=0;
            leavingModel_=0;
            injectedInModel_=0;
            changedProzessor_=0;
            noOutflowRebound_ =0;
            noRandomOutflowRebound_ = 0;
            periodicModel_ =0;
        }

        inline void reportCounters(label& particles);

        void countChangedProzessor() {changedProzessor_++;}
        void countLeaving() { leavingModel_++; }
        void countWall() { wallCollisions_++; }
        
        void countWedge() { wedgeCollisions_++; }
        void countSymmetryPlane() { symmetryPlaneCollisions_++; }

        void countInject() { injectedInModel_++; }
        void countRebound() { noOutflowRebound_++; }
        void countNoRandomOutflow() {noRandomOutflowRebound_++; }
        void countPeriodic() { periodicModel_++; }

        // Get Access to the counters
        inline label wedgeCount() const {return wedgeCollisions_;}

        const interpolation<vector>& UInterp() const {return UInterp_(); }
        const interpolation<scalar>& DInterp() const {return DInterp_(); }
        const interpolation<scalar>& rhoInterp() const {return rhoInterp_(); }
        const interpolation<vector>& gradRhoInterp() const {return gradRhoInterp_(); }
    };

protected:

    // Protected data

        //- Allow name variable access -- used for particle sampling
        nameVariableTable table_;

        // Parcel properties

        //- Supercell number
        label sCell_;
    
        //- Particle mass [kg]
        scalar mass_;
        
        //- Mass source term from the dispersed phase [kg]
        scalar dpMsource_;
    
        //- Particle weight
        scalar wt_;
    
        //- Particle drift 
        vector A_;
        vector dxDeterministic_;
        
        //- Particle random walk
        vector rw_;
        
        //- New particle position
        point toPos_;

        //- Particle Co number
        scalar maxCo_ = 0.9;


public:

    // Static data members

        //- Runtime type information
        TypeName("ItoPopeParticle");

    // Constructors

    //- Construct from owner, position, and cloud owner
    //  Other properties initialised as null
    inline ItoPopeParticle
    (
        const polyMesh& mesh,
        const barycentric& coordinates,
        const label cellI,
        const label tetFaceI,
        const label tetPtI
    );

    //- Construct from components
    inline ItoPopeParticle
    ( 
        const polyMesh& mesh,
        const barycentric& coordinates,
        const label celli,
        const label tetFacei,
        const label tetPti,
        const scalar& mass,
        const scalar& dpMsource,
        const scalar& wt,
        const vector& A
    );

    //- Construct from Istream reading field values if required
    ItoPopeParticle
    (
        const polyMesh& mesh,
        Istream& is,
        bool readFields = true,
        bool newFormat = true
    );


    //- Construct as a copy
    ItoPopeParticle(const ItoPopeParticle& p);

    //- Construct as a copy
    ItoPopeParticle(const ItoPopeParticle& p, const polyMesh& mesh);

    //- Construct and return a (basic particle) clone
    virtual autoPtr<particle> clone() const
    {
        return autoPtr<particle>(new ItoPopeParticle(*this));
    }

    //- Construct and return a (basic particle) clone
    virtual autoPtr<particle> clone(const polyMesh& mesh) const
    {
        return autoPtr<particle>(new ItoPopeParticle(*this, mesh));
    }

    //- Factory class to read-construct particles used for parallel transfer
    class iNew
    {
        const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ItoPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<ItoPopeParticle<ParticleType> >
                (
                    new ItoPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
    };

    // Member Functions

        // Access
        label& sCell()      {return sCell_;};
        label sCell() const {return sCell_;};
        scalar& m()         {return mass_;};
        scalar m() const    {return mass_;};
        scalar& dpMsource() {return dpMsource_;};
        scalar dpMsource() const        {return dpMsource_;};
        scalar& wt()        {return wt_;};
        scalar wt() const   {return wt_;};
        vector& A()         {return A_;};
        vector A() const    {return A_;};
        vector& dxDeterministic()       {return dxDeterministic_;};
        vector dxDeterministic() const  {return dxDeterministic_;};
        vector& rw()        {return rw_;};
        vector rw() const   {return rw_;};
        point& toPos()      {return toPos_;};
        point toPos() const {return toPos_;};
 

        //- Access to the name-variable table 
        nameVariableTable& nameVariableLookUpTable() {return table_;};

        const nameVariableTable& nameVariableLookUpTable() const {return table_;};

        //- Initialize settings for statistical data sampling
        void initStatisticalSampling();


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

        // Tracking

        //- Move the parcel
        template<class TrackCloudType>
        bool move(TrackCloudType& cloud, trackingData& td, const scalar trackTime);


        // Patch interactions

        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
        template<class TrackCloudType>
        bool hitPatch
        (
            TrackCloudType& cloud,
            trackingData& td,
            const label patchI,
            const scalar trackFraction,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a polyPatch
        template<class TrackCloudType>
        bool hitPatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackCloudType>
        void hitCyclicPatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a
        //  processorPatch
        template<class TrackCloudType>
        void hitProcessorPatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackCloudType>
        void hitWallPatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wedgePatch
        template<class TrackCloudType>
        void hitWedgePatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );


        //- Overridable function to handle the particle hitting a SymmetryPlanePatch
        template<class TrackCloudType>
        void hitSymmetryPlanePatch
        (
            TrackCloudType& cloud,
            trackingData& td
        );

        // Rebound from wall or noOutflow patch
        point reboundParticle
        (
            const polyPatch& pp,
            const point& toPos,
            const point& pos
        );

        //- No random outflow patch

        point noRandomOutflowRebound
        (
            const polyPatch& pp,
            const point& toPos,
            const point& pos,
            const vector& dxDet
        ) const;

        //- Periodic boundary conditions (keyword "cyclic")

        point periodicParticle
        (
            const point& toPos,
            const point& pos,
            const point& posOldBound
        ) const;

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);

        //- Set the static properties 
        template<class CloudType>
        static void setStaticProperties(CloudType& c); 

    // For Particle Sampling

        //- Return the data in a container
        //  Optionally provide a list of variables to sample
        DynamicList<scalar> getStatisticalData(const wordList& vars) const;

        //- Return variables names stored in the container
        //  Optionally provide a list of variables to sample
        DynamicList<word> getStatisticalDataNames(const wordList& vars) const;
        

        // I-O

        //- Read particle data for all particles in a cloud from files
        template<class CloudType>
        static void readFields(CloudType& c);

        //- Write particle data for all particles in a cloud to a file
        template<class CloudType>
        static void writeFields(const CloudType& c);


        // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const ItoPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ItoPopeParticleI.H"
#include "ItoPopeParticleTrackingDataI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ItoPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
