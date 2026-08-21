/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ThermoPopeCloud

Description
    Templated base class for thermodynamic cloud
    Creates the thermodynamic model used for the chemistryModel
    of the particles

    - Adds to Ito cloud
    - Heat transfer

    Includes 
    - ReactionThermo model
    - SLGThermo 

SourceFiles
    ThermoPopeCloudI.H
    ThermoPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoPopeCloud_H
#define ThermoPopeCloud_H

#include "ItoPopeCloud.H"
#include "thermoPopeCloud.H"

#include "SLGThermo.H"
#include "fvMatrices.H"

#include "mmcVarSet.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes

template<class CloudType>
class CompositionModel;

template<class CloudType>
class ThermoPhysicalCouplingModel;

/*---------------------------------------------------------------------------*\
                    Class ThermoPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class ReactionThermo>
class ThermoPopeCloud
:
    public CloudType,
    public thermoPopeCloud
{
    public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of the reaction thermo model
        typedef ReactionThermo ReactionThermoType;

        //- Type of particle the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef ThermoPopeCloud<CloudType,ReactionThermo> ThermoPopeCloudType;

    private:

    // Private data

        //- Cloud copy pointer
        autoPtr<ThermoPopeCloud<CloudType,ReactionThermo> > cloudCopyPtr_;

        //- Reaction thermo model used for particle chemsitry model
        ReactionThermo* thermo_{nullptr};

        //- Boolean to check if the thermo model is found in the objectRegistry
        //  or was created by the ThermoPopeCloud
        bool createdThermoModel_{false};

    // Private member functions

        //- Disallow default bitwise copy construct
        ThermoPopeCloud(const ThermoPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ThermoPopeCloud&);


    protected:

    // Protected data

        //- SLG thermodynamics package
        autoPtr<SLGThermo> SLGThermo_;

        //- Temperature [K]
        const volScalarField& T_;

        //- Pressure [Pa]
        const volScalarField& p_;
    

    // Sub Models 

        //- Composition sub-model
        autoPtr<CompositionModel<ThermoPopeCloudType > >  compositionModel_;
        
        //- Density coupling model 
        autoPtr<ThermoPhysicalCouplingModel<ThermoPopeCloudType > > thermoPhysicalCouplingModel_;

    // Protected Member Functions

        //- Check that size of a composition field is valid
        void checkSuppliedComposition
        (
            const scalarField& YSupplied,
            const scalarField& Y,
            const word& YName
        );

    // Initialisation

        //- Set cloud sub-models
        void setModels(const mmcVarSet& Xi);

        //- Reset state of cloud
        void cloudReset(ThermoPopeCloudType& c);

        //- Lookup or construct thermo model if not found
        ReactionThermo* lookupOrConstructThermo();

    public:

    // Constructors

        //- Construct with given fields
        ThermoPopeCloud
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
        ThermoPopeCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const mmcVarSet& Xi,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

    //- Destructor
    virtual ~ThermoPopeCloud();


    // Member Functions

    // Access

        //- Return a reference to the cloud copy
        inline const ThermoPopeCloud& cloudCopy() const;

        //- Return const access to thermo package
        inline const ReactionThermo& thermo() const;
        
        //- Return non-const access to thermo package
        inline ReactionThermo& thermo();
        
        //- Return const access to SLG thermo package
        inline const SLGThermo& slgThermo() const;

        //- Return const access to the carrier temperature field
        inline const volScalarField& T() const;

        //- Return const access to the carrier prressure field
        inline const volScalarField& p() const;

    // Sub-models

        //- Return const access to composition model
        inline const CompositionModel<ThermoPopeCloudType>& composition() const;

        //- Return const access to density coupling model
        inline const ThermoPhysicalCouplingModel<ThermoPopeCloudType>& coupling() const;
        

	//- Second-conditioning query: always false at this level.
        //  MixingPopeCloud overrides this with the real check.
        virtual bool secondCondMixingEnabled() const { return false; }

        //- Return non-const access to composition model
        inline CompositionModel<ThermoPopeCloudType>& composition() ;

        //- Return non-const access to density coupling model
        inline ThermoPhysicalCouplingModel<ThermoPopeCloudType>& coupling();

    // Sources

    // Check

        //- Maximum temperature
        inline scalar Tmax() const;

        //- Minimum temperature
        inline scalar Tmin() const;

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
        
    // I-O

        virtual void writeFields() const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

