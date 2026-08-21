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
    Foam:mixParticleModel

Description
    Abstract base class for particle interaction models
    
    - builds the particle list and searches the closest pairs
    - mixing of the two particles is then described by the chosen mixing model

    Settings are controlled over the mixingModel coeff dictionary. E.g. for 
    the MMCCurl model:
    MMCCurlCoeffs
    {
        ri          1.0E-3;     // Normalisation factor for the reference space
        
        Ximi                    // Reference variable normalisation
        {
          fm                       0.03;
        }

        pairingMethod              subVolumes;  // pairing method
                                                // Options: local,global, subVolumes

        // Optional: Exclude particles from parallel mixing
        // Particles with a reference variable less than fLow or higher than
        // fHigh are only mixed locally on the processor
        fLow             0.0001;
        fHigh             0.99;
    }

    Details of the mixing models can be found in: 
    [1] Son Vo, O.T. Stein, A. Kronenburg, M.J. Cleary,
        Assessment of mixing time scales for a sparse particle method,
        Combustion and Flame,
        Volume 179,
        2017,
        Pages 280-299,
        ISSN 0010-2180,
        https://doi.org/10.1016/j.combustflame.2017.02.017. 

\*---------------------------------------------------------------------------*/

#ifndef mixParticleModel_H
#define mixParticleModel_H

#include "CloudMixingModel.H"
#include <vector>
#include <list>
//#include "kdTree.H"
#include "particlePtrList.H"
#include "particlePairingMethod.H"
#include "baseParticleDataContainer.H"
#include "mixingSubVolumes.H"
//#include "particleMatchingAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mixParticleModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class mixParticleModel
:
    public CloudMixingModel<CloudType>
{
    typedef typename CloudType::particleType particleType;
    typedef Tuple2<scalarList,particleType> mixParticle;
    
    public:
    
        // Struct to store the eulerian field properties for each particle
        // Required for the aISO model 
        class eulerianFieldData
        : public baseParticleDataContainer
        {
            public:
                // Constructors 
                eulerianFieldData()
                : baseParticleDataContainer(2,7,2) {};
                
                //- Construct from input stream
                eulerianFieldData(Istream& is)
                : baseParticleDataContainer(is) 
                {};
                
            // Non-Const Access
                //- random number
                inline scalar& Rand() {return scalarValues_[0];}
                
                //- Effective diffusivity 
                inline scalar& DEff() {return scalarValues_[1];}
                
                //- Laminar diffusivity
                inline scalar& D() {return scalarValues_[2];}
                
                //- Turbulent diffusivity
                inline scalar& Dt() {return scalarValues_[3];}
                
                //- Cell size of LES mesh
                inline scalar& DeltaE() {return scalarValues_[4];}
                
                //- Viscosity
                inline scalar& mu() {return scalarValues_[5];}
                
				inline scalar& vb() {return scalarValues_[6];}

                //- Index of the particleList_
                inline label& particleIndex() {return labelValues_[0];}
        
                inline label& processorIndex() {return labelValues_[1];}
                
                inline List<scalar>& magSqrRefVar() {return scalarFields_[0];}
                
                inline scalarField& XiR() {return scalarFields_[1];}

            // Const Access
                inline const scalar& Rand() const {return scalarValues_[0];}
                inline const scalar& DEff() const {return scalarValues_[1];}
                inline const scalar& D() const {return scalarValues_[2];}
                inline const scalar& Dt() const {return scalarValues_[3];}
                inline const scalar& DeltaE() const {return scalarValues_[4];}
                inline const scalar& mu() const {return scalarValues_[5];}
				inline const scalar& vb() const {return scalarValues_[6];}
                inline const label& particleIndex() const {return labelValues_[0];}
                inline const label& processorIndex() const {return labelValues_[1];}
                inline const List<scalar>& magSqrRefVar() const {return scalarFields_[0];}
                inline const scalarField& XiR() const {return scalarFields_[1];}
        };
        
        // Required for kd-Tree like construction of finding particle pairs
        class lessArg
        {
            label ii_=-1;
            
            public:
        
                lessArg(label ii) {ii_=ii;} 

                inline bool operator()  
                (       
                    const eulerianFieldData& one,
                    const eulerianFieldData& two
                ) const
                {
                    if (ii_ < 3)
                    {
                        return one.position()[ii_] < two.position()[ii_]; 
                    }
                    else 
                    {
                        label k = ii_-3;
                        return one.XiR()[k]< two.XiR()[k];
                    }

                    // Return a default to avoid warnings of the compiler
                    FatalError << "Comparison operator did not work properly"<<exit(FatalError);
                    return false;
                };
        };
            
    protected:
        
        // Reference Variable Settings
            
            //- Information about the reference variable 
            const mmcVarSet& XiR_;

            //- Names of the reference variables
            const wordList XiRNames_;
            
            //- Number of reference variables
            const label numXiR_;      

            //- Normalisation factors
	    //  Default size matches numXiR_ (one per first-conditioning
            //  reference variable). Derived models that use a different
            //  reference space (e.g. secondCondMMCcurl with 4D
            //  (phiModified, xi_x, xi_y, xi_z)) overwrite this list in
            //  their constructor.
            List<scalar> Xii_;
            const scalar ri_;
            //const List<scalar> Xii_;

        // (Optional) Discard particles for parallel mixing if their reference 
        // variable, XiR(), is below fLow or above fMax

            //- Lower limit of the reference variable (default: 0.0)
            //const scalar fLow_;

            //- Upper limit of the reference variable (default: 1.0)
            //const scalar fHigh_;

           
        // Const References to Eulerian volScalarFields

            //- Reference to effective diffusivity
            const volScalarField& DEff_;
        
            //- Reference to laminar diffusivity
            const volScalarField& D_;

            //- Reference to turbulent diffusivity
            const volScalarField& Dt_;
            
            //- Eulerian grid size
            const volScalarField& DeltaE_;

        // field needed for soot 
            const volScalarField& mu_;

			const volScalarField& vb_;

        // required to get processors for mixing

            //- Stores the method selected for particle pairing
            const particlePairingMethod pairingMethod_;

            mixingSubVolumes mixSubVolumes_;

            //autoPtr<particleMatchingAlgorithm<eulerianFieldData>>
            //    particlePairAlgorithm_;

        // Internal data

            //- List of all particles stored as pointers
            particlePtrList<particleType> particleList_;
        
            //- List of eulerian data fields for each particle 
            DynamicList<eulerianFieldData> eulerianFieldDataList_;
            
            //- List of particle pairs
            //  Pair indices access eulerianFieldDataList_ the index of the 
            //  particle can be found with particleIndex() of the 
            //  eulerianField
            DynamicList<List<label>> particlePairs_;

            //- Store a copy of all particles from neigh processors
            PtrDynList<particleType> particleListProcs_;
            
            //- Contains processor IDs of all processors in which particles 
            //- will mix
            List<label> particleMixingProcessors_;

            //- List of starting indices for the particle data of each processor
            List<label> startIndexOfParticle_;


	    //- Diagnostic counter: number of times each XiR axis was selected
            //  as the split axis in KkdTreeLikeSearch during the most recent
            //  findPairs() call. Indexed 0..numXiR_-1 (per process).
            //  Reset at the start of findPairs(). Mutable so it can be written
            //  from the const k-d tree routine.
            mutable List<label> splitAxisHistogram_;

        // Private Member Functions
        
            //- Build particle list and correct for parallel runs
            //void buildParticleList();

	    //  Virtual so derived models can filter the particle set before pairing
            virtual void buildParticleList();

            //- Build particle list for particles considered only for local 
            //  mixing
            //  See also fLow_ and fHigh_
            //  --> Note: This clears the old particleList_ 
            //void buildParticleListLocalMixing(const scalar fLow, const scalar fHigh);

            //- Get the processor IDs in which particles will mix
            List<label> getParticleMixingProcessors();

            //- Correct particle list for parallel runs 
            void correctParticleListParallel();

            //- Retrieve Eulerian data fields from processors
            //  Required for mixing and finding the pairing
            void collectEulerianDataFields();

            //- Get the particle data fields for mixing
            void collectParticleData();

            //- Mix the list of particles
            virtual void SmixList() ;

            //- Read from the dictionary the normalisation factors for the
            //- reference variable
            List<scalar> getXiNormalisation();
            
            void findPairs
            (
                const DynamicList<eulerianFieldData>& eulerianFieldList,
                DynamicList<List<label>>& pairs
            ) const;
            
            void KkdTreeLikeSearch
            (
                const DynamicList<eulerianFieldData>& particleList,
                label l,
                label u,
                std::vector<label>& pInd,
                std::vector<label>& L,
                std::vector<label>& U
            ) const;

    public:
        // Constructors

            //- Construct from dictionary
            mixParticleModel
            (
                const dictionary& dict,
                CloudType& owner,
                const word& type,
                const mmcVarSet& Xi
            );

            //- Construct copy
            mixParticleModel 
            (
                const mixParticleModel<CloudType>& cm
            );
            

        // Destructor
        
            //- Default destructor
            virtual ~mixParticleModel() = default;


        // Member functions

            //- General mixing function
            virtual void Smix();
               
            
            //- Function mixing particle pairs
            //  Avoid direct call - use Smix()
            virtual void mixpair
            (
                particleType& p,
                const eulerianFieldData& pEulerianFields,
                particleType& q,
                const eulerianFieldData& qEulerianFields,
                scalar& deltaT
            ) = 0;
               
               
        // Inline
               
            // inside mmcVarSet
            inline label numXiR() const
            {
                return numXiR_;
            };

	    //- Per-process histogram of split-axis selections from the most
            //  recent findPairs() invocation. Index i corresponds to XiR[i].
            inline const List<label>& splitAxisHistogram() const
            {
                return splitAxisHistogram_;
            };

            //- Read-only access to the assembled particle list and the pair
            //  list (useful for derived-class diagnostics).
            inline const particlePtrList<particleType>& particleList() const
            {
                return particleList_;
            };

            inline const DynamicList<List<label>>& particlePairs() const
            {
                return particlePairs_;
            };

            // inside mmcVarSet
            inline const wordList& XiRNames() const
            {
                return XiRNames_;
            };
            
            inline const mmcVarSet& XiR() const
            {
                return XiR_;
            };
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixParticleModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

