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
    Foam::MixingPopeCloud

Description
    Templated base class for mixing cloud

SourceFiles
    MixingPopeCloudI.H
    MixingPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef MixingPopeCloud_H
#define MixingPopeCloud_H

#include "ItoPopeCloud.H"
#include "ThermoPopeCloud.H"
#include "mixingPopeCloud.H"
#include "SLGThermo.H"
#include "OUProcessUpdate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class CloudMixingModel;

/*---------------------------------------------------------------------------*\
                    Class MixingPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MixingPopeCloud
:
    public CloudType,
    public mixingPopeCloud
{
    public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of particle the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef MixingPopeCloud<CloudType> mixingPopeCloudType;

    private:

    // Private data

        //- Cloud copy pointer
        autoPtr<MixingPopeCloud<CloudType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        MixingPopeCloud(const MixingPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const MixingPopeCloud&);

    protected:

    // Protected data

        //- Mixing model
	//- Primary (first-conditioning) mixing model
        autoPtr<CloudMixingModel<MixingPopeCloud<CloudType> > >  mixingModel_;

	//- Second-conditioning mixing model (optional; null when disabled)
        autoPtr<CloudMixingModel<MixingPopeCloud<CloudType> > >  secondCondMixingModel_;

	// Second conditioning parameters (read from cloudProperties)
 
            //- Fraction of particles assigned to the second-cond subset [0,1]
            scalar secondCondR_;
 
            //- OU weight β used in φ° = φ·exp(β·ω_OU)
            scalar secondCondBeta_;
 
            //- OU process relaxation time τ_OU [s]
            scalar secondCondTauOU_;

	    //- Unburnt mixture temperature Tu [K]
            scalar secondCondTu_;
 
            //- Burnt mixture temperature Tb [K]
            scalar secondCondTb_;
 
            //- Pre-exponential coefficient in the progress-variable reaction term
            scalar secondCondAPhi_;
 
            //- Activation parameter in the progress-variable reaction term
            scalar secondCondZPhi_;

    // Protected Member Functions

        //- Set cloud sub-models
        void setModels(const mmcVarSet& Xi);

        //- Reset state of cloud
        void cloudReset(MixingPopeCloud<CloudType>& c);


  public:

    // Constructors

        //- Construct with given fields
        MixingPopeCloud
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
        MixingPopeCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const mmcVarSet& Xi,
            const Switch initAtCnstr = true,
            bool readFields = true        
        );

    //- Destructor
    virtual ~MixingPopeCloud();


    // Member Functions

    // Access

        //- Return a reference to the cloud copy
        inline const MixingPopeCloud& cloudCopy() const;


    // Sub-models

        //- Return const access to mixing model
        inline const CloudMixingModel<MixingPopeCloud<CloudType> >& mixing() const;

        //- Return non-const access to mixing model
        inline CloudMixingModel<MixingPopeCloud<CloudType> >& mixing();

	//- True if the second-conditioning mixing model has been constructed
        inline bool secondCondMixingEnabled() const;
 
        //- Return non-const access to second-conditioning mixing model
        //  (only valid when secondCondMixingEnabled() == true)
        inline CloudMixingModel<MixingPopeCloud<CloudType> >& secondCondMixing();
 
        //- Return const access to second-conditioning mixing model
        inline const CloudMixingModel<MixingPopeCloud<CloudType> >& secondCondMixing() const;

    // Second conditioning
 
        //- Fraction of particles in the second-conditioning subset
        inline scalar secondCondR() const;
 
        //- OU weight β
        inline scalar secondCondBeta() const;
 
        //- OU relaxation time τ_OU [s]
        inline scalar secondCondTauOU() const;

	//- Unburnt temperature Tu [K]
        inline scalar secondCondTu() const;
 
        //- Burnt temperature Tb [K]
        inline scalar secondCondTb() const;
 
        //- Pre-exponential reaction coefficient A_phi
        inline scalar secondCondAPhi() const;
 
        //- Activation parameter Z_phi
        inline scalar secondCondZPhi() const;

    // Second-conditioning cloud-level operations
    // (called from moveParticles.H in the correct timestep sequence)
 
        //- Apply W(φ) = A·(1−φ)·exp[Z·(φ−1)] to ALL particles for one Δt.
        //  Must be called after mixing().Smix() and before updateOUProcess().
        void updatePhiReaction(scalar deltaT);
 
        //- Advance ω_OU and recompute φ° = φ·exp(β·ω_OU) for flagged particles.
        //  Must be called after updatePhiReaction() and before secondCondMixing().Smix().
        void updateOUProcess(scalar deltaT);

    // I-O
    
        //- Add new particle

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
};


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MixingPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MixingPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

