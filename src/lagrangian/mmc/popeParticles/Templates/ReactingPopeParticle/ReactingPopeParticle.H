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
    Foam::ReactingPopeParticle

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes Kinematic parcel sub-models, plus:
    - heat transfer

SourceFiles
    ReactingPopeParticleI.H
    ReactingPopeParticle.C
    ReactingPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingPopeParticle_H
#define ReactingPopeParticle_H

#include "particle.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class ReactingPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const ReactingPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class ReactingPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ReactingPopeParticle
:
    public ParticleType
{


public:
    
    typedef ReactingPopeParticle<ParticleType> reactingParticleIOType;
    
    class trackingData
    :
        public ParticleType::trackingData
    {
        private:

        // Private data

        public:

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                TrackCloudType& cloud
            )
            :
                ParticleType::trackingData(cloud)
            {}
    };

private:
    
    // Private member variables
    
        //- CPU time required to solve the reaction 
        scalar cpuTime_{0};

public:

    // Static data members

        //- Runtime type information
        TypeName("ReactingPopeParticle");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ReactingPopeParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        ReactingPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        ReactingPopeParticle(const ReactingPopeParticle& p);

        //- Construct as a copy
        ReactingPopeParticle(const ReactingPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ReactingPopeParticle(*this, mesh));
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

            autoPtr<ReactingPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<ReactingPopeParticle<ParticleType> >
                (
                    new ReactingPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt,
                const label cellI
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt,
                const label cellI
            );
            
            //- Mix operator for this particleType
            static void mixProperties
            (
                ReactingPopeParticle<ParticleType>& p, 
                ReactingPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            static void mixProperties
            (
                ReactingPopeParticle<ParticleType>& p, 
                ReactingPopeParticle<ParticleType>& q, 
                const scalar& mixExtent,
                const scalar& mixExtentSoot
            );

            //- Mix operator for this particleType
            static void mixProperties
            (
                ReactingPopeParticle<ParticleType>& p, 
                ReactingPopeParticle<ParticleType>& q, 
                scalar mixExtent,
                scalarList ScaledExtent
            );

        // Access to member variables
        
            //- CPU time
            inline scalar& cpuTime();
            
            //- const access to CPU time
            inline const scalar& cpuTime() const;


        // I-O

            //- Read particle data for all particles in a cloud from files
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write particle data for all particles in a cloud from files
            template<class CloudType>
            static void writeFields(const CloudType& c);

    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const ReactingPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingPopeParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
