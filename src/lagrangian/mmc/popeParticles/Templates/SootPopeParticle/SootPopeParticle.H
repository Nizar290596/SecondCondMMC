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
    Foam::SootPopeParticle

Description


SourceFiles
    SootPopeParticleI.H
    SootPopeParticle.C
    SootPopeParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef SootPopeParticle_H
#define SootPopeParticle_H

#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class SootPopeParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const SootPopeParticle<ParticleType>&
);

/*---------------------------------------------------------------------------*\
                       Class SootPopeParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class SootPopeParticle
:
    public ParticleType
{
public:
    
    typedef SootPopeParticle<ParticleType> sootParticleIOType;
    
    class trackingData
    :
        public ParticleType::trackingData
    {
        public:

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                TrackCloudType& cloud
            );
    };

protected:

    // Protected data

            //- Soot volume fraction
            scalarField fVsoot_;

            //- Intermittency []
            scalarField Intermittency_;

            //- Soot volme fraction [dimless] 
            scalar sootVf_;

            //- Soot intermittency [dimless] 
            scalar sootIm_;

            //- soot mass fraction [dimless] 
            scalar ySoot_;

            //- soot number density [1/kg-mixture] 
            scalar nSoot_;

            //- soot contributions of each step 
            //  it will be optimized later, currently a list of scalar is used
            //  -> soot mass fraction:
            //     ySootInC2H2, ySootSgC2H2, ySootOxO2, ySootOxOH
            //  -> soot number density:
            //     nSootInC2H2, nSootCg
            scalarList sootContribs_;

            //- Soot contribution names
            static const wordList sootContribsNames_;

            //- XO2/XN2 (used for exp. comparison)
            scalar XO2XN2_;


    // Protected Member Functions

public:

    // Static data members

        //- Runtime type information
        TypeName("SootPopeParticle");

        //- String representation of properties
        AddToPropertyList
        (
            ParticleType,
            " fvSoot"
          + " Intermittency"
          + " sootVf" 
          + " sootIm" 
          + " ySoot"  
          + " nSoot"  
          + " sootContribs" 
          + " XO2XN2"       
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline SootPopeParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from Istream
        SootPopeParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        SootPopeParticle(const SootPopeParticle& p);

        //- Construct as a copy
        SootPopeParticle(const SootPopeParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new SootPopeParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new SootPopeParticle(*this, mesh));
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

            autoPtr<SootPopeParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<SootPopeParticle<ParticleType> >
                (
                    new SootPopeParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to soot volume fraction
            inline const scalarField& fVsoot() const;

            //- Return const access to Intermittency
            inline const scalarField& Intermittency() const;

            //- Return const access to soot volume fraction 
            inline const scalar& sootVf() const;

            //- Return const access to soot intermittency 
            inline const scalar& sootIm() const;

            //- Return const access to soot mass fraction 
            inline const scalar& ySoot() const;

            //- Return const access to soot number density 
            inline const scalar& nSoot() const;

            //- Return const access to soot conributions from each steps 
            inline const scalarList& sootContribs() const;

            //- Return const access to XO2XN2 
            inline const scalar& XO2XN2() const;


        // Edit

            //- Return access to soot volume fraction
            inline scalarField& fVsoot();

            //- Return const access to Intermittency
            inline scalarField& Intermittency();


            //- Return access to soot volume fraction 
            inline scalar& sootVf();

            //- Return access to soot intermittency 
            inline scalar& sootIm();
            
            //- Return access to soot mass fraction 
            inline scalar& ySoot();

            //- Return access to soot number density [1/kg-mixture] 
            inline scalar& nSoot();

            //- Return access to soot conributions from each steps 
            inline scalarList& sootContribs();

            //- Return access to XO2XN2 
            inline scalar& XO2XN2();

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
            //  This mixing does not consider soot modeling
            static void mixProperties
            (
                SootPopeParticle<ParticleType>& p, 
                SootPopeParticle<ParticleType>& q, 
                scalar mixExtent
            );

            //- Mix operator for this particleType
            //  Considers soot model, see sootMMCCurl model
            static void mixProperties
            (
                SootPopeParticle<ParticleType>& p, 
                SootPopeParticle<ParticleType>& q, 
                const scalar& mixExtent,
                const scalar& mixExtentSoot
            );


            static void mixProperties
            (
                SootPopeParticle<ParticleType>& p, 
                SootPopeParticle<ParticleType>& q, 
                scalar mixExtent,
                scalarList ScaledExtent
            );
        
        // For particle sampling 

            //- Initialize settings for statistical data sampling
            void initStatisticalSampling();

            //- Return the data in a container
            DynamicList<scalar> getStatisticalData
            (
                const wordList& vars
            ) const;

            //- Return variables names stored in the container
            DynamicList<word> getStatisticalDataNames
            (
                const wordList& vars
            ) const;

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
            const SootPopeParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SootPopeParticleI.H"
#include "SootPopeParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SootPopeParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
