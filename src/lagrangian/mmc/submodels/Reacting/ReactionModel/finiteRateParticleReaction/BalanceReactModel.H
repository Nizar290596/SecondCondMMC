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
    Foam::BalanceReactModel

Description
    A class for balancing load on reaction

\*---------------------------------------------------------------------------*/

#ifndef BalanceReactModel_H
#define BalanceReactModel_H

#include "IOmanip.H"
#include "ReactionModel.H"
#include "baseParticleDataContainer.H"
#include <vector>
#include <list>
#include <unordered_map> 
#include <chrono>

//#include "ParticleTuple.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class BalanceReactModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BalanceReactModel
:
    public ReactionModel<CloudType>
{
    private:
    
        // particle data container for storing and sending reduced particle 
        // information
        class particleDataContainer
        : public baseParticleDataContainer
        {            
            public:
                
                particleDataContainer()
                : baseParticleDataContainer(2,5,1) {}
            
                particleDataContainer(Istream& is)
                : baseParticleDataContainer(is) 
                {};
            
                // non const access
                inline scalar& p() {return scalarValues_[0];}
                inline scalar& T() {return scalarValues_[1];}
                inline scalar& h() {return scalarValues_[2];}
                inline scalar& mixtureFraction() {return scalarValues_[3];}
                inline scalar& cpuTime() {return scalarValues_[4];}
                
                inline label& proc() {return labelValues_[0];}
                inline label& origId() {return labelValues_[1];}
                
                inline scalarField& Y() {return scalarFields_[0];}
            
                // const access
                inline scalar p() const {return scalarValues_[0];}
                inline scalar T() const {return scalarValues_[1];}
                inline scalar h() const {return scalarValues_[2];}
                inline scalar mixtureFraction() const {return scalarValues_[3];}
                inline scalar cpuTime() const {return scalarValues_[4];}
                
                inline label proc() const {return labelValues_[0];}
                inline label origId() const {return labelValues_[1];}
                
                inline const scalarField& Y() const {return scalarFields_[0];}
        };

    protected:

        //- Switch to check if it is called the first time in the simulation
        bool firstTime_{true};

        scalar totalCpuTime_;
        
    // Protected Member functions
        
        //- Build the particleDataContainer list from the cloud particles
        DynamicList<particleDataContainer> buildReactPaticleList();
        
        //- solve the reaction for all particles in the given list
        void solveReaction
        (
            UList<particleDataContainer>& parList
        );

        //- Calculate the particles to send/recv to/from other 
        //- processors
        //  Returns a list of sublists of particles to send to processor i
        void particlesToSend
        (
            const DynamicList<particleDataContainer>& parList,
            const scalar cpuTimeToSend,
            const label& start,
            label& end
        );
        
        //- Get percentage of processor data to send to other processors
        //- to balance the processor load
        Tuple2
        <
            List<Tuple2<scalar,label>>,
            List<label>
        > 
        getProcessorBalancing();

        //- Get cpu times of each processor in a descending order
        //  Format:
        //  List[index].first:  CPU time required for last time step
        //  List[index].second: Returns an Foam::Pair<label> with
        //                        Pair.first():  processor ID
        //                        Pair.second(): Number of particles on that processor 
        List<std::pair<scalar,Pair<label>>> getSortedCPUTimesOnProcessor() const;
        
        //- Update the totalCpuTime_ variable 
        void updateTotalCpuTime
        (
            const DynamicList<particleDataContainer>& reactParList
        );

    public:


    // Constructors

        //- Construct from dictionary
        BalanceReactModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

    // Destructor
        virtual ~BalanceReactModel();

    // Public Member Functions 
    
        //- Balanced solving of the reaction 
        void SreactBalance();
        
        //- Calculate chemistry
        //  Here only a dummy function
        virtual void calculate
        (
            const scalar t0, 
            const scalar dt,
            const scalar ha, 
            const scalar pc, 
            const scalar pT,
            const scalar pz,
            scalarField& Y
        ) override
        {Info << "dummy calculate function"<<endl;}
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BalanceReactModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
