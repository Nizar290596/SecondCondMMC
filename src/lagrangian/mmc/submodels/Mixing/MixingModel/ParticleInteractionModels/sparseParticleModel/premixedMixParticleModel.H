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
    Foam:premixedMixParticleModel

Description
    Abstract base class for particle interaction models
    
    - builds the particle list and searches the closest pairs
    - mixing of the two particles is then described by the chosen mixing model

    Settings are controlled over the mixingModel coeff dictionary. E.g. for 
    the premixedMMCCurl model:
    premixedMMCCurlCoeffs
    {
        ...
        ri          1.0E-3;     // Normalisation factor for the reference space
        mixAllProcs true;       // Consider all processors for mixing
                                // default is true
        maxRefDist  0.5;        // Optional paramter only needed if mixAllProcs
                                // is set to false!
                                // Meaning: particles with a larger value than
                                // maxRefDist might be more suitable mixing 
                                // candidates if they exist
    }

    Details of the mixing models can be found in: 
    [1] Iaroslavtceva, N., Kronenburg, A. & Stein, O.T,
        Multiple Mapping Conditioning Mixing Time Scales for Turbulent Premixed Flames,
        Flow Turbulence Combust,
        Volume 110,
        2023,
        Pages 395–415,
        https://doi.org/10.1007/s10494-022-00375-1. 

\*---------------------------------------------------------------------------*/

#ifndef premixedMixParticleModel_H
#define premixedMixParticleModel_H

#include "CloudMixingModel.H"
#include <vector>
#include <list>
#include "kdTree.H"
#include "particlePtrList.H"
#include "particlePairingMethod.H"
#include "baseParticleDataContainer.H"
#include "mixingSubVolumes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class premixedMixParticleModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class premixedMixParticleModel
:
    public CloudMixingModel<CloudType>
{
    typedef typename CloudType::particleType particleType;
    
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
                
                //- Thickening factor from ATF
                inline scalar& TF() {return scalarValues_[6];}
                
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
                inline const scalar& TF() const {return scalarValues_[6];}
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
            const scalar ri_;
            const List<scalar> Xii_;

           
        // Const References to Eulerian volScalarFields
        
            //- Reference to laminar diffusivity
            const volScalarField& D_;

            //- Reference to turbulent diffusivity
            const volScalarField& Dt_;
            
            //- Eulerian grid size
            const volScalarField& DeltaE_;

           // field needed for soot 
            const volScalarField& mu_;
            
        // Required when ATF is used
            
            // thickening factor field needed for kdTree search
            const volScalarField& TF_;
            
            // true: flame thickening is acounted for in kdTree pairing (cLESm = const/TF)
            // false: standard kdTree pairing (cLESm = const)
            const Switch ATFkd_;

        // required to get processors for mixing

            //- Stores the method selected for particle pairing
            const particlePairingMethod pairingMethod_;

            mixingSubVolumes mixSubVolumes_;

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

        // Private Member Functions
        
            //- Build particle list
            void buildParticleList();

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
            
            //- Find pairs of particles
            void findPairs
            (
                const DynamicList<eulerianFieldData>& eulerianFieldList,
                DynamicList<List<label>>& pairs
            ) const;
            
            //- Search for particle pairs by using a recursive call to a 
            //- kdTree construction with box size 2
            void premixedkdTreeLikeSearch
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
            premixedMixParticleModel
            (
                const dictionary& dict,
                CloudType& owner,
                const word& type,
                const mmcVarSet& Xi
            );

            //- Construct copy
            premixedMixParticleModel 
            (
                const premixedMixParticleModel<CloudType>& cm
            );
            

        // Destructor
        
            //- Default destructor
            virtual ~premixedMixParticleModel() = default;


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
#   include "premixedMixParticleModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
