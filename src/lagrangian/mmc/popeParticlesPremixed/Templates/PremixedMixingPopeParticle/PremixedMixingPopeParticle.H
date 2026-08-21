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
    Foam::PremixedMixingPopeParticle

Description


SourceFiles
    PremixedMixingPopeParticleI.H
    PremixedMixingPopeParticle.C
    PremixedMixingPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef PremixedMixingPopeParticle_H
#define PremixedMixingPopeParticle_H

#include "particle.H"
#include "mmcVarSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class PremixedMixingPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const PremixedMixingPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class PremixedMixingPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class PremixedMixingPopeParticle
:
    public ParticleType
{
public:
    
    typedef PremixedMixingPopeParticle<ParticleType> mixingParticleIOType;
    
    class trackingData
    :
        public ParticleType::trackingData
    {
        private:

        // Private data


            // Interpolators for continuous phase fields
            
                //- Reference variable list of field interpolators
                PtrList<interpolation<scalar> > XiRInterp_;

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
            //  reference variables
            inline const PtrList<interpolation<scalar> >& XiRInterp() const;
    };
    
protected:

    // Protected data
    
    // Static data
        static HashTable<label, word> indexInXiR_;
        
        //- Name of the reference variables
        static List<word> XiRNames_;

        // Particle properties

            //- Reference mixture fraction
            scalarField XiR_;

            //- Distance between mixing particle pairs in reference space
            scalarField dXiR_;
            
            //- Distance between mixing particle pairs in physical space
            scalar dx_;
            
            //- Mixing time scale of a particle
            scalar mixTime_;
            
            //- Mixing extent of a particle pair
            scalar mixExt_;
            
            //- Progress variable interpolated from the Eulerian field 
            scalar progVarInt_;
            
            //- Progress variable relaxation factor
            scalar XiRrlxFactor_;
            
            //- Progress variable relaxation step (counting restarts when the relaxation factor changes)
            scalar rlxStep_;
            
    // Protected Member Functions


public:

    // Static data members

        //- Runtime type information
        TypeName("PremixedMixingPopeParticle");

        //- String representation of properties
        AddToPropertyList
        (
            ParticleType,
            " XiR"
          + " dXiR"
          + " dx"
          + " mixTime"
          + " mixExt"
          + " progVarInt"
          + " XiRrlxFactor"
          + " rlxStep"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline PremixedMixingPopeParticle
        (
            const polyMesh& mesh,
            //const vector& position,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        PremixedMixingPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        PremixedMixingPopeParticle(const PremixedMixingPopeParticle& p);

        //- Construct as a copy
        PremixedMixingPopeParticle(const PremixedMixingPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new PremixedMixingPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new PremixedMixingPopeParticle(*this, mesh));
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

            autoPtr<PremixedMixingPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<PremixedMixingPopeParticle<ParticleType> >
                (
                    new PremixedMixingPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

	    //- Second conditioning is not used in the premixed path;
            //  return 1 so all particles are treated as chemistry-active.
            inline label secondCondFlag() const { return 1; }

            
            //- Return const access to reference variables
            inline const scalarField& XiR() const;

            //- Return const access to the distance between mixing particles in Xi-space
            inline const scalarField& dXiR() const;

            //- Return const access to the distance between mixing particles in x-space
            inline scalar dx() const;
            
            //- Return const access to specific reference variable
            inline const scalar& XiR(const word& name) const;
            
            //- Return const access to specific ref variable mixing distance 
            inline const scalar& dXiR(const word& name) const;
            
            //- Return const access to the mixing time scale of a particle
            inline scalar mixTime() const;
            
            //- Return const access to the mixing extent of a particle pair
            inline scalar mixExt() const;
            
            //- Return const access to the progress variable interpolated from the Eulerian field
            inline scalar progVarInt() const;
            
            //- Return const access to the progress variable relaxation factor
            inline scalar XiRrlxFactor() const;
            
            //- Return const access to the progress variable relaxation step
            inline scalar rlxStep() const;


        // Edit
            
            //- Return reference variables
            inline scalarField& XiR();

            //- Return distance between mixing particles in Xi-space
            inline scalarField& dXiR();
            
            //- Return distance between mixing particles in x-space
            inline scalar& dx();
            
            //- Return specific reference variable
            inline scalar& XiR(const word& name);
            
            //- Return specific reference variable mixing distance
            inline scalar& dXiR(const word& name);
            
            //- Return mixing time scale of a particle
            inline scalar& mixTime();
            
            //- Return mixing extent of a particle pair
            inline scalar& mixExt();
            
            //- Return progress variable interpolated from the Eulerian field
            inline scalar& progVarInt();
            
            //- Return progress variable relaxation factor
            inline scalar& XiRrlxFactor();
            
            //- Return progress variable relaxation step
            inline scalar& rlxStep();
            
            //- Set the static member variables 
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
            
            // Mix operator for this particleType
            static void mixProperties
            (
                PremixedMixingPopeParticle<ParticleType>& p, 
                PremixedMixingPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            static void mixProperties
            (
                PremixedMixingPopeParticle<ParticleType>& p, 
                PremixedMixingPopeParticle<ParticleType>& q, 
                const scalar& mixExtent,
                const scalar& mixExtentSoot
            );

            static void mixProperties
            (
                PremixedMixingPopeParticle<ParticleType>& p, 
                PremixedMixingPopeParticle<ParticleType>& q, 
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

            //- Read
            template<class CloudType, class MixingType>
            static void readFields
            (
                CloudType& c,
                const MixingType& mixModel //Needed to access reference variables
            );
            
            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class MixingType>
            static void writeFields
            (
                const CloudType& c,
                const MixingType& mixModel //Needed to access reference variables
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);

            
    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const PremixedMixingPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PremixedMixingPopeParticleI.H"
#include "PremixedMixingPopeParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PremixedMixingPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
