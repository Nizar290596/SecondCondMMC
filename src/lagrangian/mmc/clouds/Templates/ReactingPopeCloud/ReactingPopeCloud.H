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
    Foam::ReactingPopeCloud

Description
    Templated base class for thermodynamic cloud

    - Adds to kinematic cloud
    - Heat transfer

SourceFiles
    ReactingPopeCloudI.H
    ReactingPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingPopeCloud_H
#define ReactingPopeCloud_H

#include "ItoPopeCloud.H"
#include "reactingPopeCloud.H"
#include "fluidThermo.H"
#include "mmcVarSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//template<class CloudType>
//class HeatTransferModel;

template<class CloudType>
class ReactionModel;

/*---------------------------------------------------------------------------*\
Class ReactingPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingPopeCloud
:
    public CloudType,
    public reactingPopeCloud
{
    public:

        // Public typedefs

            //- Type of cloud this cloud was instantiated for
            typedef CloudType cloudType;

            //- Type of parcel the cloud was instantiated for
            typedef typename CloudType::particleType particleType;

            //- Convenience typedef for this cloud type
            typedef ReactingPopeCloud<CloudType> reactingPopeCloudType;


    private:

        // Private data

            //- Cloud copy pointer
            autoPtr<ReactingPopeCloud<CloudType> > cloudCopyPtr_;

            // Private member functions

            //- Disallow default bitwise copy construct
            ReactingPopeCloud(const ReactingPopeCloud&);

            //- Disallow default bitwise assignment
            void operator=(const ReactingPopeCloud&);

            //-swritch of balance reaction load
            Switch balanceReactionLoad_; 


    protected:

        // Protected data
        
        // References to the cloud sub-models

            //- Reaction model
            autoPtr<ReactionModel<ReactingPopeCloud<CloudType> > > reactionModel_;

        // Sources

        // Protected Member Functions

        // Initialisation

        //- Set cloud sub-models
        void setModels();

        // Cloud evolution functions

        //- Reset state of cloud
        void cloudReset(ReactingPopeCloud<CloudType>& c);


    public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingPopeCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& DEff,
            const volScalarField& rho,
            const volVectorField& gradRho,
            const mmcVarSet& Xi,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

        //- Construct and lookup required fields in the object registry
        ReactingPopeCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const mmcVarSet& Xi,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

        //- Destructor
        virtual ~ReactingPopeCloud();

        // Member Functions

        // Access

        // Sub-models

            //- Return const access to reacting model
            inline const ReactionModel<ReactingPopeCloud<CloudType> >& reaction() const;

            //- Return const access to reacting model
            inline ReactionModel<ReactingPopeCloud<CloudType> >& reaction();

            //- Return const access to balanceReactionLoad swtich
            inline const Switch& balanceReactionLoad();

        // Modelling options

        // Sources

        // Enthalpy

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

        //- Solve the cloud - calls all evolution functions
        template<class TrackCloudType> 
        void solve
        (
            TrackCloudType& cloud,
            typename particleType::trackingData& td
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
