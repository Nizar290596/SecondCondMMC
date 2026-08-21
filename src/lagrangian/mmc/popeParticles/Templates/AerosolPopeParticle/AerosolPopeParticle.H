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
    Foam::AerosolPopeParticle

Description
    particle which introduces the aerosol method

SourceFiles
    AerosolPopeParticleI.H
    AerosolPopeParticle.C
    AerosolPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef AerosolPopeParticle_H
#define AerosolPopeParticle_H

#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CloudType>
class AerosolPopeCloud;

template<class ParticleType>
class AerosolPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const AerosolPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class AerosolPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class AerosolPopeParticle
:
public ParticleType
{
    // Friend declarations
        template<class CloudType> 
        friend class AerosolPopeCloud;  

    // Static data
            
            static List<word> physicalAerosolPropertyNames_;
            
            static List<word> psdPropertyNames_;

protected:    
    // Protected data
    
        // Particle properties

            //- physical properties of the aerosol 
            scalarField physicalAerosolProperties_;
    
            //- PSD properties
            scalarField psdProperties_;
 
    // Protected Member Functions
                        
public:
    
    typedef AerosolPopeParticle<ParticleType> aerosolParticleIOType;
    
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
            );
    };
    
    // Static data members

        //- Runtime type information
        TypeName("AerosolPopeParticle");

        //- String representation of properties
        AddToPropertyList
        (
            ParticleType,
          + " physicalAerosolProperties"
          + " psdProperties"
        );

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline AerosolPopeParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        AerosolPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        AerosolPopeParticle(const AerosolPopeParticle& p);

        //- Construct as a copy
        AerosolPopeParticle(const AerosolPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new AerosolPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new AerosolPopeParticle(*this, mesh));
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

            autoPtr<AerosolPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<AerosolPopeParticle<ParticleType> >
                (
                    new AerosolPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions
        // Access
        
            //- Return const access to physical aerosol properties
            inline const scalarField& physicalAerosolProperties() const;
            
            //- Return const access to PSD properties used for the aerosol method 
            inline const scalarField& psdProperties() const;

        // Edit
            
            //- Return access to physical aerosol properties
            inline scalarField& physicalAerosolProperties();

            //- Return access to PSD properties used for the aerosol method
            inline scalarField& psdProperties();

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
            
            // Mix operator for this particleType
            static void mixProperties
            (
                AerosolPopeParticle<ParticleType>& p, 
                AerosolPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            static void mixProperties
            (
                AerosolPopeParticle<ParticleType>& p, 
                AerosolPopeParticle<ParticleType>& q, 
                const scalar& mixExtent,
                const scalar& mixExtentSoot
            );

            static void mixProperties
            (
                AerosolPopeParticle<ParticleType>& p, 
                AerosolPopeParticle<ParticleType>& q, 
                scalar mixExtent,
                scalarList ScaledExtent
            );
            
            //- Set static particle properties
            template<class CloudType>
            static void setStaticProperties(CloudType& c); 

        // For particle sampling 

            //- Initialize settings for statistical data sampling
            void initStatisticalSampling();

            //- Return the data in a container
            DynamicList<scalar> getStatisticalData(const wordList& vars) const;

            //- Return variables names stored in the container
            DynamicList<word> getStatisticalDataNames(const wordList& vars) const;

        // I-O

            //- Read
            template<class CloudType, class AerosolType>
            static void readFields
            (
                CloudType& c,
                const AerosolType& aerosolModel    
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class AerosolType>
            static void writeFields
            (
                const CloudType& c,
                const AerosolType& aerosolModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);

    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const AerosolPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AerosolPopeParticleI.H"
#include "AerosolPopeParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AerosolPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
