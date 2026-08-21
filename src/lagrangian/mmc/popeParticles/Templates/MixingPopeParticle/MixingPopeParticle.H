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
    Foam::MixingPopeParticle

Description


SourceFiles
    MixingPopeParticleI.H
    MixingPopeParticle.C
    MixingPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef MixingPopeParticle_H
#define MixingPopeParticle_H

#include "particle.H"
#include "mmcVarSet.H"
#include "OUProcessUpdate.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CloudType>
class MixingPopeCloud;

template<class ParticleType>
class MixingPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const MixingPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class MixingPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class MixingPopeParticle
:
    public ParticleType
{
    // Friend declarations
        template<class CloudType> 
        friend class MixingPopeCloud;  

public:
    
    typedef MixingPopeParticle<ParticleType> mixingParticleIOType;
    
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

	//- Cached β from cloud (used in mixProperties where no cloud ref available)
        //static scalar secondCondBeta_s_;

        // Particle properties

            //- Reference mixture fraction
            scalarField XiR_;

            //- Distance between mixing particle pairs in reference space
            scalarField dXiR_;
            
            //- Distance between mixing particle pairs in physical space
            scalar dx_;

	   // Second conditioning fields
 
                //- Subset flag: 0 = regular particle, 1 = second-cond subset
                label secondCondFlag_;
 
                //- Ornstein-Uhlenbeck process state variable ω_OU
                scalar omegaOU_;

		//- Progress variable φ = (T - Tu)/(Tb - Tu)  ∈ [0,1]
                scalar phi_;
 
                //- Modified progress variable φ° = φ * exp(β * ω_OU)
                scalar phiModified_;
public:

    // Static data members

        //- Runtime type information
        TypeName("MixingPopeParticle");

        //- String representation of properties
        AddToPropertyList
        (
            ParticleType,
            " XiR"
          + " dXiR"
          + " dx"
	  + " secondCondFlag"
	  + " omegaOU"
	  + " phi"
	  + " phiModified"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline MixingPopeParticle
        (
            const polyMesh& mesh,
            //const vector& position,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        MixingPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        MixingPopeParticle(const MixingPopeParticle& p);

        //- Construct as a copy
        MixingPopeParticle(const MixingPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new MixingPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new MixingPopeParticle(*this, mesh));
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

            autoPtr<MixingPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<MixingPopeParticle<ParticleType> >
                (
                    new MixingPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access
            
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

	    //- Return const access to second conditioning flag
            inline label secondCondFlag() const;
 
            //- Return const access to OU state variable
            inline scalar omegaOU() const;

	    //- Return const access to progress variable
            inline scalar phi() const;
 
            //- Return const access to modified progress variable
            inline scalar phiModified() const;

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

	    //- Return reference to second conditioning flag
            inline label& secondCondFlag();
 
            //- Return reference to OU state variable
            inline scalar& omegaOU();

	    //- Return reference to progress variable
            inline scalar& phi();
 
            //- Return reference to modified progress variable
            inline scalar& phiModified();
            
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
                MixingPopeParticle<ParticleType>& p, 
                MixingPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            static void mixProperties
            (
                MixingPopeParticle<ParticleType>& p, 
                MixingPopeParticle<ParticleType>& q, 
                const scalar& mixExtent,
                const scalar& mixExtentSoot
            );

            static void mixProperties
            (
                MixingPopeParticle<ParticleType>& p, 
                MixingPopeParticle<ParticleType>& q, 
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
            const MixingPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MixingPopeParticleI.H"
#include "MixingPopeParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MixingPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

