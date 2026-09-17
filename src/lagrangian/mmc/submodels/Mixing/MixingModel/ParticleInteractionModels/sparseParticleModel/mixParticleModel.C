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
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY s
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "mixParticleModel.H"
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "processorPolyPatch.H"
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::mixParticleModel<CloudType>::mixParticleModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const mmcVarSet& Xi
)
:
    CloudMixingModel<CloudType>(dict,owner,type),
  
    XiR_(Xi),

    XiRNames_(Xi.rVarInXi().toc()),
    
    numXiR_(XiRNames_.size()),
    
    ri_(readScalar(this->coeffDict().lookup("r_i"))),

    Xii_(getXiNormalisation()),
    physicalLocalization_(this->coeffDict().lookupOrDefault("physicalLocalization",
        owner.cloudProperties().subOrEmptyDict("secondConditioning").template lookupOrDefault<bool>("enabled", false))),
    retainedPairFraction_(this->coeffDict().lookupOrDefault("retainedPairFraction",
        owner.cloudProperties().subOrEmptyDict("secondConditioning").template lookupOrDefault<bool>("enabled", false) ? scalar(0.8) : scalar(1))),
    maxPairDistance_(this->coeffDict().lookupOrDefault("maxPairDistance", scalar(0))),
    mixingExtentModel_(this->coeffDict().template lookupOrDefault<word>("mixingExtentModel",
        type == "secondCondMMCcurl" ? word("modifiedCurl") : word("exponential"))),
    mixingTimeScale_(this->coeffDict().template lookupOrDefault<word>("mixingTimeScale",
        this->coeffDict().lookupOrDefault("aISO", true) ? word("aISO") : word("gradient"))),
    prescribedTauMix_(this->coeffDict().lookupOrDefault("tauMix", scalar(-1))),

//    fLow_(this->coeffDict().template lookupOrDefault<scalar>("fLow",-GREAT)),
    
//    fHigh_(this->coeffDict().template lookupOrDefault<scalar>("fHigh",GREAT)),
    
    DEff_(owner.mesh().objectRegistry::lookupObject<volScalarField>("DEff")),

    D_(owner.mesh().objectRegistry::lookupObject<volScalarField>("D")),

    Dt_(owner.mesh().objectRegistry::lookupObject<volScalarField>("Dt")),

    DeltaE_
    (
        owner.mesh().objectRegistry::lookupObject<volScalarField>("DeltaE")
    ),

    mu_
    (
        owner.mesh().objectRegistry::lookupObject<volScalarField>("thermo:mu")
    ),

    vb_(owner.mesh().objectRegistry::lookupObject<volScalarField>("vb")),
    
    pairingMethod_(this->coeffDict()),

    mixSubVolumes_(owner.mesh(),this->coeffDict(),ri_,Xii_[0])

    //particlePairAlgorithm_
    //(
    //    particleMatchingAlgorithm<eulerianFieldData>::New
    //    (
    //        ri_,
    //        Xii_,
    //        this->coeffDict()
    //    )
    //)
{
    if (!(ri_ > 0 && std::isfinite(ri_) && retainedPairFraction_ > 0
        && retainedPairFraction_ <= 1 && maxPairDistance_ >= 0
        && std::isfinite(maxPairDistance_)))
        FatalErrorInFunction << "Require r_i > 0, 0 < retainedPairFraction <= 1, "
            << "and finite maxPairDistance >= 0" << exit(FatalError);
    if (mixingExtentModel_ != "exponential" && mixingExtentModel_ != "modifiedCurl")
        FatalErrorInFunction << "mixingExtentModel must be exponential or modifiedCurl" << exit(FatalError);
    if (mixingTimeScale_ != "aISO" && mixingTimeScale_ != "prescribed" && mixingTimeScale_ != "gradient")
        FatalErrorInFunction << "mixingTimeScale must be aISO, gradient or prescribed" << exit(FatalError);
    if (mixingTimeScale_ == "prescribed" && !(prescribedTauMix_ > 0 && std::isfinite(prescribedTauMix_)))
        FatalErrorInFunction << "Prescribed mixing requires finite tauMix > 0 [s]" << exit(FatalError);
    for (const word key : {word("nPairSamples"), word("particleFilter"),
                           word("localnessLimited"), word("fullSort")})
        if (this->coeffDict().found(key))
            WarningInFunction << key << " is not used; remove it. Use the explicit "
                << "localization and retainedPairFraction controls." << nl;
    Info<< type << ": physicalLocalization=" << physicalLocalization_
        << ", retainedPairFraction=" << retainedPairFraction_
        << ", maxPairDistance=" << maxPairDistance_
        << ", extent=" << mixingExtentModel_ << ", timescale=" << mixingTimeScale_
        << ", tauMix=" << prescribedTauMix_ << nl;
}


template <class CloudType>
Foam::mixParticleModel<CloudType>::mixParticleModel
(
    const mixParticleModel<CloudType>& cm
)
:
    CloudMixingModel<CloudType>(cm),
    
    XiR_(cm.XiR_),

    XiRNames_(cm.XiRNames_),
    
    numXiR_(XiRNames_.size()),

    ri_(readScalar(this->coeffDict().lookup("r_i"))),

    Xii_(cm.Xii_),
    physicalLocalization_(cm.physicalLocalization_),
    retainedPairFraction_(cm.retainedPairFraction_),
    maxPairDistance_(cm.maxPairDistance_),
    mixingExtentModel_(cm.mixingExtentModel_),
    mixingTimeScale_(cm.mixingTimeScale_),
    prescribedTauMix_(cm.prescribedTauMix_),

//    fLow_(cm.fLow_),
    
//    fHigh_(cm.fHigh_),

    DEff_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>("DEff")
    ),  

    //initialise the fields needed for aISO
    D_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>("D")
    ),

    Dt_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>("Dt")
    ),

    DeltaE_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>
        ("DeltaE")
    ),

    mu_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>("thermo:mu")
    ),
    
    vb_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("vb")),
    

    pairingMethod_(cm.pairingMethod_),
    
    mixSubVolumes_(cm.mixSubVolumes_)

    //particlePairAlgorithm_(cm.particlePairAlgorithm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::mixParticleModel<CloudType>::buildParticleList
(
    //const scalar fLow,
    //const scalar fHigh
)
{
    // ========================================================================
    // Build local particleList

    PtrList<volScalarField> magSqr_XiR_(XiRNames_.size());
    
    label II = 0;
    
    for (const word& nameI : XiRNames_)
    {
		//Info << "XIRNAMES_ " << XiRNames_ <<endl;
		//Info << "nameI: "<<nameI << endl;
        magSqr_XiR_.set
        (
            this->XiR_.rVarInXiR()[nameI],
            new volScalarField
            (
                magSqr(fvc::grad(this->XiR_.Vars(nameI).field()))
            )
        );
        II++;
    }
    
    PtrList<interpolationCellPoint<scalar> > magSqr_intp_(magSqr_XiR_.size());
    
    forAll(magSqr_XiR_, vsfI)
    {   
        magSqr_intp_.set
        (
            vsfI,
            new interpolationCellPoint<scalar>(magSqr_XiR_[vsfI])
        );
    }
            
    interpolationCellPoint<scalar> DEff_intp_(this->DEff_);

    //for aISO
    interpolationCellPoint<scalar> D_intp_(this->D_);

    interpolationCellPoint<scalar> Dt_intp_(this->Dt_);

    interpolationCellPoint<scalar> DeltaE_intp_(this->DeltaE_);

    interpolationCellPoint<scalar> mu_intp_(this->mu_);

    interpolationCellPoint<scalar> vb_intp_(this->vb_);
    
  
    // clear particle list from old data
    particleList_.clear();
    
    StochasticLib1& rand = this->owner().rndGen();
    
    eulerianFieldDataList_.clear();
    
    // running index for particle position
    label particleInd=0;
    
    forAllIters(this->owner(), iter)
    {
        // Only add particle to the list if it is in the range of 
        // flow to fMax
        //if (iter().XiC().first() < fLow || iter().XiC().first() > fHigh)
        //   continue;

        // Asign pointer to particle to list of particles
        particleList_.append
        (
            iter.get()
        );
        
        // Store all eulerian fields 
        eulerianFieldData eulerianFields;
        
        // Get the position, cell and face of the particle
        const vector pos = iter().position();
        const label cellI = iter().cell();
        const label faceI = iter().face();
        
        eulerianFields.particleIndex() = particleInd++;
        
        eulerianFields.processorIndex() = Pstream::myProcNo();
        
        // The eulerian data field also has to store the position for the 
        // k-d tree later
        eulerianFields.position() = pos;
        
        // Also store the reference variable
        eulerianFields.XiR() = iter().XiR();
        
        eulerianFields.Rand() = rand.Random();
        
        eulerianFields.DEff() = DEff_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.D() = D_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.Dt() = Dt_intp_.interpolate(pos,cellI,faceI);
            
        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos,cellI,faceI);
            
        eulerianFields.mu() = mu_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.vb() = vb_intp_.interpolate(pos,cellI,faceI);
        eulerianFields.magSqrRefVar().resize(iter().XiR().size());
        // Reference Variables & related quantitites 
        forAll(iter().XiR(), j)
        {
			eulerianFields.XiR()[j] = iter().XiR()[j];
            eulerianFields.magSqrRefVar()[j] = magSqr_intp_[j].interpolate
            (
                pos,cellI,faceI
            );            
        }
            
        eulerianFieldDataList_.append
        (
            std::move(eulerianFields)
        );
        
    }
    
	
    //if (eulerianFieldDataList_.size() < 2)
    //{
	//	Info<<"EulerianField Data List Cleared!"<<endl;
    //    particlePairs_.clear();
     //   return;
    //}


    // if run in parallel get all required particles of neighbouring processors
    if 
    (
            Pstream::parRun() 
         && pairingMethod_.method() != particlePairingMethod::localPairing
    )
    {
        // findPairs is called in correctParticleListParallel
        //Info << "Start CorrectParticleListParallel " << endl;
        correctParticleListParallel();
    }
    else
    {
        findPairs(eulerianFieldDataList_,particlePairs_);   
    }
}


/*template <class CloudType>
void Foam::mixParticleModel<CloudType>::buildParticleListLocalMixing
(
    const scalar fLow,
    const scalar fHigh
)
{
    // ========================================================================
    // Build local particleList

    PtrList<volScalarField> magSqr_XiR_(XiRNames_.size());
    
    label II = 0;
    
    for (const word& nameI : XiRNames_)
    {
        magSqr_XiR_.set
        (
            this->XiR_.rVarInXiR()[nameI],
            new volScalarField
            (
                magSqr(fvc::grad(this->XiR_.Vars(nameI).field()))
            )
        );
        II++;
    }
    
    PtrList<interpolationCellPoint<scalar> > magSqr_intp_(magSqr_XiR_.size());
    
    forAll(magSqr_XiR_, vsfI)
    {   
        magSqr_intp_.set
        (
            vsfI,
            new interpolationCellPoint<scalar>(magSqr_XiR_[vsfI])
        );
    }
            
    interpolationCellPoint<scalar> DEff_intp_(this->DEff_);

    //for aISO
    interpolationCellPoint<scalar> D_intp_(this->D_);

    interpolationCellPoint<scalar> Dt_intp_(this->Dt_);

    interpolationCellPoint<scalar> DeltaE_intp_(this->DeltaE_);

    interpolationCellPoint<scalar> mu_intp_(this->mu_);

    interpolationCellPoint<scalar> vb_intp_(this->vb_);
  
    // clear particle list from old data
    particleList_.clear();
    
    StochasticLib1& rand = this->owner().rndGen();
    
    eulerianFieldDataList_.clear();
    
    // running index for particle position
    label particleInd=0;
    
    forAllIters(this->owner(), iter)
    {
        // Only add particle to the list if they are not already considered 
        // in the parallel mixing
        if (iter().XiC().first() >= fLow && iter().XiC().first() <= fHigh)
            continue;

        // Asign pointer to particle to list of particles
        particleList_.append
        (
            iter.get()
        );
        
        // Store all eulerian fields 
        eulerianFieldData eulerianFields;
        
        // Get the position, cell and face of the particle
        const vector pos = iter().position();
        const label cellI = iter().cell();
        const label faceI = iter().face();
        
        eulerianFields.particleIndex() = particleInd++;
        
        eulerianFields.processorIndex() = Pstream::myProcNo();
        Info << "eulFeildsProcIndx" << eulerianFields.processorIndex() << endl;
        
        // The eulerian data field also has to store the position for the 
        // k-d tree later
        eulerianFields.position() = pos;
        
        // Also store the reference variable
        eulerianFields.XiR() = iter().XiR();
        
        eulerianFields.Rand() = rand.Random();
        
        eulerianFields.DEff() = DEff_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.D() = D_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.Dt() = Dt_intp_.interpolate(pos,cellI,faceI);
            
        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos,cellI,faceI);
            
        eulerianFields.mu() = mu_intp_.interpolate(pos,cellI,faceI);
        eulerianFields.vb() = vb_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.magSqrRefVar().resize(iter().XiR().size());
        // Reference Variables & related quantitites 
        forAll(iter().XiR(), j)
        {
            eulerianFields.magSqrRefVar()[j] = magSqr_intp_[j].interpolate
            (
                pos,cellI,faceI
            );            
        }
            
        eulerianFieldDataList_.append
        (
            std::move(eulerianFields)
        );
    }
    
    if (eulerianFieldDataList_.size() < 2)
    {
        particlePairs_.clear();
        return;
    }

    particlePairAlgorithm_->findPairs(eulerianFieldDataList_,particlePairs_);
}*/




template<class CloudType>
void Foam::mixParticleModel<CloudType>::correctParticleListParallel()
{  
    particleMixingProcessors_ = getParticleMixingProcessors();
    //Info << "particleMixingProcessors " << particleMixingProcessors_ << endl;

    // Collect the eulerian data from other processors
    //Info << "Start CollectEulerianDataFields " << endl;
    collectEulerianDataFields();
//	Info << "End CorrectParticleListParallel " << endl;
    // Find the particle pairs to mix
    //Info << "Start FindingPairs " << endl;
    findPairs(eulerianFieldDataList_,particlePairs_);
    //Info << "End FindingPairs " << endl;  
    collectParticleData();

    // Important note: The particle pointers cannot all be stored in 
    // particleList_ because this is a special class that does not 
    // deallocate objects when clear() is called. The PtrDynList however
    // does deallocate once clear() is called. 
    // Also this has to be outside the loop for the processors, as 
    // reallocation of the dynamic list would invalidate the references 
    // to the particle locations
    forAll(particleListProcs_,i)
    {
        particleList_.append(&particleListProcs_[i]);
    }
    //Info << "ParticlePairs" << particlePairs_ << endl;
}


template<class CloudType>
void Foam::mixParticleModel<CloudType>::collectEulerianDataFields()
{
    // First send all the eulerian data fields and create the pairing lists
    // Then only send pairs that mix with particles located on the current
    // processor 
    //Info << "CheckPoint 1 " <<endl;
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
	//Info << "CheckPoint 2 " << endl;
    // send the list of particles of this processor to neighbour 
    // processors
    for (auto& procI : particleMixingProcessors_)
    {
        if (procI != Pstream::myProcNo())
        {
            UOPstream toBuffer(procI,pBufs);
            toBuffer << eulerianFieldDataList_.size();
            //Info<<"List content " << eulerianFieldDataList_ << endl;
            for (auto& e : eulerianFieldDataList_)
                toBuffer << e;
        }
    }
    //Info << "CheckPoint 3 " <<endl;
    pBufs.finishedSends();
    //Info << "CheckPoint 4 " <<endl;
    startIndexOfParticle_.clear();
    startIndexOfParticle_.resize(Pstream::nProcs(),-1);
    label previousParticleSize = 0;
//	Info << "CheckPoint 5 " <<endl;
    // Loop over all processors to update the startIndexOfParticle list
    // it is important that all processors have the same order of Eulerian
    // fields. Otherwise they might calculate different pairings!
//	Info << "CheckPoint 4 " <<endl;
    DynamicList<eulerianFieldData> tlocalEulerianFields = 
        std::move(eulerianFieldDataList_);
    eulerianFieldDataList_.clear();
//	Info << "CheckPoint 5 " <<endl;
    // Estimate space for eulerianFields
    eulerianFieldDataList_.reserve
    (
        particleMixingProcessors_.size()*tlocalEulerianFields.size()
    );
//	Info << "CheckPoint 6 " <<endl;
    for (const label& procI : particleMixingProcessors_)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream fromBuffer(procI,pBufs);
            label size;
            fromBuffer >> size;

            eulerianFieldDataList_.reserve
                (eulerianFieldDataList_.size()+size);
            
            startIndexOfParticle_[procI] = previousParticleSize;
            previousParticleSize += size;

            for (label k=0; k < size; k++)
                eulerianFieldDataList_.append
                (
                    eulerianFieldData(fromBuffer)
                );
        }
        else
        {
            startIndexOfParticle_[procI] = previousParticleSize;
            previousParticleSize += tlocalEulerianFields.size();
            eulerianFieldDataList_.append(std::move(tlocalEulerianFields));
        }
    }
    //Info << "CheckPoint 7 " <<endl;
}


template<class CloudType>
void Foam::mixParticleModel<CloudType>::collectParticleData()
{
    // ====================================================================
    //         Note to sending particles between processors 
    // ====================================================================
    // The particle class cannot be moved or has a copy assignment
    // constructor, as the field const fvMesh& mesh_; cannot be moved
    // or copied. The only way to transfer particles between processors
    // is to follow the approach of Cloud::move() where particles are 
    // streamed to the buffer and then read from it using the constructor:
    // particleType(mesh,is);
    // ====================================================================

    // Go over the pairs and check which particles mix with particles 
    // located on the current processor.
    // Keep track of the particles that need to be send to other processors
    List<DynamicList<label>> particlesToSendToProcessor(Pstream::nProcs());


    // Const reference to the current mesh
    const fvMesh& mesh = this->owner().mesh();

    // particle number on local processor
    const label numLocalParticles = particleList_.size();

    // Reserve some space
    for (label i : particleMixingProcessors_)
        particlesToSendToProcessor[i].reserve(0.05*numLocalParticles);
    
    for (const List<label>& pair : particlePairs_)
    {
        // Does the pair contain a particle from another processor
        bool processorParticle = false;
        // Does the pair contain a local particle
        bool localParticle = false;
        
        DynamicList<label> localParticleToSend(pair.size());
        DynamicList<label> processorToSend(pair.size());


        for (const label& i : pair)
        {
            if (eulerianFieldDataList_[i].local())
            {
                localParticleToSend.append(i);
                localParticle = true;
            }
            else
            {
                processorToSend.append
                (
                    eulerianFieldDataList_[i].processorIndex()    
                );
                processorParticle = true;
            }
        }
        
        if (localParticle && processorParticle)
        {
            for (label i : localParticleToSend)
            {
                for (label& procI : processorToSend)
                {
                    particlesToSendToProcessor[procI].append(i);
                }
            }
        }
    }


    for (label proc : particleMixingProcessors_)
    {
        auto& send = particlesToSendToProcessor[proc];
        std::sort(send.begin(), send.end());
        const auto end = std::unique(send.begin(), send.end());
        send.resize(std::distance(send.begin(), end));
    }
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (const label& procI : particleMixingProcessors_)
    {
        if (procI != Pstream::myProcNo())
        {
            UOPstream toBuffer(procI,pBufs);
            
            toBuffer << particlesToSendToProcessor[procI].size()<<nl;

            // Stream all particles and their eulerian field to the buffer
            forAll(particlesToSendToProcessor[procI],i)
            {
                const eulerianFieldData& e1 = 
                    eulerianFieldDataList_[particlesToSendToProcessor[procI][i]];
                toBuffer << e1.particleIndex();
                toBuffer << particleList_[e1.particleIndex()];
            }
        }
    }

    pBufs.finishedSends();
    
    // Store the particles in the particleListProcs_ and add a reference 
    // to them to the particleList_
    particleListProcs_.clear();
    
    label particleInd = numLocalParticles;

    for (const label& procI : particleMixingProcessors_)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream fromBuffer(procI,pBufs);

            label nParticlesToRead = 0;
            fromBuffer >> nParticlesToRead;

            // reserve space
            particleListProcs_.reserve
            (
                particleListProcs_.size()+nParticlesToRead
            );

            // Read all particles and eulerianFields
            for (label k=0; k < nParticlesToRead; k++)
            {
                label procParticleInd;
                fromBuffer >> procParticleInd;

                particleListProcs_.append
                (
                    new particleType(mesh,fromBuffer)
                );

                eulerianFieldDataList_
                [
                    startIndexOfParticle_[procI]+procParticleInd
                ].particleIndex() = particleInd++;
            }
        }
    }
}


template<class CloudType>
void Foam::mixParticleModel<CloudType>::Smix()
{
    // First mix particles considered for local mixing
    //buildParticleListLocalMixing(fLow_,fHigh_);
    // Mix the particles
    //SmixList();
    // Now all particles for which parallel handling is considered
    minMixingTime_ = GREAT;
    maxMixingCourant_ = 0;
    weightedProgressVarianceLoss_ = 0;
    buildParticleList();
    SmixList();
    const label interval = this->coeffDict().template lookupOrDefault<label>("diagnosticInterval", 100);
    if (this->owner().mesh().time().timeIndex() % max(interval, label(1)) == 0)
    {
        scalar minimum = minMixingTime_, maximum = maxMixingCourant_;
        reduce(minimum, minOp<scalar>());
        reduce(maximum, maxOp<scalar>());
        scalar loss = weightedProgressVarianceLoss_, weight = 0;
        forAllIters(this->owner(), particle) weight += particle().wt();
        reduce(loss, sumOp<scalar>());
        reduce(weight, sumOp<scalar>());
        Info<< "Mixing timescale: min tau=" << (minimum < GREAT ? minimum : 0)
            << ", max exchange dt/tau=" << maximum
            << ", extent=" << mixingExtentModel_ << nl;
        if (loss > 0 && weight > VSMALL)
            Info<< "Dense progress mixing: Nphi="
                << loss/(2*this->owner().mesh().time().deltaTValue()*weight)
                << " [1/s], inferred from weighted variance loss during mixing" << nl;
    }
}


template<class CloudType>
Foam::List<Foam::label>
Foam::mixParticleModel<CloudType>::getParticleMixingProcessors()
{
    if (pairingMethod_.global())
    {
        List<label> procList(Pstream::nProcs());
        forAll(procList,i)
        {
            procList[i] = i;
        }
        return procList;
    }

    if (pairingMethod_.method() == particlePairingMethod::neighbourPairs)
    {
        if (!processorGraphReady_)
        {
            List<labelList> neighbours(Pstream::nProcs());
            DynamicList<label> localNeighbours;
            const polyBoundaryMesh& patches = this->owner().mesh().boundaryMesh();
            forAll(patches, i)
                if (isA<processorPolyPatch>(patches[i]))
                    localNeighbours.append(refCast<const processorPolyPatch>(patches[i]).neighbProcNo());
            neighbours[Pstream::myProcNo()] = localNeighbours;
            Pstream::gatherList(neighbours);
            Pstream::scatterList(neighbours);
            forAll(neighbours, rank)
                for (label other : neighbours[rank])
                    if (other != rank)
                        processorEdges_.emplace_back(min(rank, other), max(rank, other));
            std::sort(processorEdges_.begin(), processorEdges_.end());
            processorEdges_.erase(std::unique(processorEdges_.begin(), processorEdges_.end()), processorEdges_.end());
            processorGraphReady_ = true;
        }
        const auto partners = secondConditioningNumerics::neighbourPartners
            (Pstream::nProcs(), processorEdges_, this->owner().mesh().time().timeIndex());
        const label me = Pstream::myProcNo(), other = partners[me];
        List<label> group(other < 0 ? 1 : 2);
        group[0] = other < 0 ? me : min(me, other);
        if (other >= 0) group[1] = max(me, other);
        return group;
    }

    // Only works for one reference variable 
    if (numXiR_ > 1)
        FatalError << "More than one reference variable selected."<<nl
            << "Only particle pairing local and global are possible"
            << exit(FatalError);

    
    // get the sub-volume
    return mixSubVolumes_.getSubVolume(Pstream::myProcNo());
}


template<class CloudType>
void Foam::mixParticleModel<CloudType>::SmixList()
{
    scalar deltaT = this->owner().mesh().time().deltaT().value();

    for (auto& pair : particlePairs_)
    {
        if(pair.size() == 2)
        {
            const eulerianFieldData& e1 = eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = eulerianFieldDataList_[pair[1]];

            // Only call mixpair for processor local particles, as remote
            // particles are mixed on their respective processor
            if (e1.local() || e2.local())
                mixpair
                (
                    particleList_[e1.particleIndex()],e1,
                    particleList_[e2.particleIndex()],e2,
                    deltaT
                );
        }
        else if(pair.size() == 3)
        {
            const eulerianFieldData& e1 = eulerianFieldDataList_[pair[0]];
            const eulerianFieldData& e2 = eulerianFieldDataList_[pair[1]];
            const eulerianFieldData& e3 = eulerianFieldDataList_[pair[2]];

            if (e1.local() || e2.local() || e3.local())
            {
                scalar halfDt = 0.5*deltaT;
                mixpair(particleList_[e1.particleIndex()],e1,
                        particleList_[e2.particleIndex()],e2,halfDt);
                mixpair(particleList_[e2.particleIndex()],e2,
                        particleList_[e3.particleIndex()],e3,halfDt);
                mixpair(particleList_[e3.particleIndex()],e3,
                        particleList_[e1.particleIndex()],e1,halfDt);
            }
        }
    }
}


template <class CloudType>
List<scalar> Foam::mixParticleModel<CloudType>::getXiNormalisation() 
{
    // dictionary to read the normalisation parameters for 
    // the reference variables 
    const dictionary XiDict(this->coeffDict().subDict("Xim_i"));

    Info << nl << "The Ximi parameters are: "<< XiDict << endl;

    List<scalar> Xii(numXiR_);

    const HashTable<label, word>& XiRIndexes = this->XiR().rVarInXiR();

    label i=0;
    for (const word& refVarName :this->XiRNames())
    {
        const label index = XiRIndexes[refVarName];
        Xii[index] = XiDict.found(refVarName+"_m")
            ? readScalar(XiDict.lookup(refVarName+"_m")) : scalar(1);
        if (!XiDict.found(refVarName+"_m")
            && this->coeffDict().lookupOrDefault("includeShadowPositions", true))
            FatalErrorInFunction << "Missing normalization " << refVarName << "_m" << exit(FatalError);
        if (!(Xii[index] > 0 && std::isfinite(Xii[index])))
            FatalErrorInFunction << "Invalid normalization for " << refVarName << exit(FatalError);
    }

	//Info << "Xii" << Xii << endl;
    return Xii;
}

template<class CloudType>
void Foam::mixParticleModel<CloudType>::findPairs
(
    const DynamicList<eulerianFieldData>& eulerianFieldList,
    DynamicList<List<label>>& pairs
) const
{
    // Clear particle pairs first
    pairs.clear();

    // Reset the per-axis split-counter (one slot per XiR axis)
    splitAxisHistogram_.setSize(Xii_.size());
    splitAxisHistogram_ = 0;

    if (eulerianFieldList.size() < 2) return;

    // Keeping track of indices for premixedkdTreeLikeSearch
    std::vector<label> L;
    std::vector<label> U;
    L.reserve(eulerianFieldList.size());
    U.reserve(eulerianFieldList.size());
    
    // create an index list for the particle data
    std::vector<label> pInd(eulerianFieldList.size());
    std::iota(pInd.begin(),pInd.end(),0);
    
    KkdTreeLikeSearch(eulerianFieldList,1,eulerianFieldList.size(),pInd,L,U);

    // reserve space for list of pairs
    pairs.reserve(ceil(0.5*eulerianFieldList.size()));

    for(size_t i=0; i<L.size(); i++)
    {
        label p = L[i] - 1;

        label q = L[i];

        if(U[i] - L[i] < 2)
        {
            List<label> pair(2);
            pair[0] = pInd[p];
            pair[1] = pInd[q];
            
            pairs.append(std::move(pair));
        }
        else if(U[i] - L[i] == 2)
        {
            label r = L[i] + 1;
            
            List<label> pair(3);
            pair[0] = pInd[p];
            pair[1] = pInd[q];
            pair[2] = pInd[r];

            pairs.append(std::move(pair));
        }
	//Info << "PAIRS SIZE" << pairs.size() << endl;
    }

    if (retainedPairFraction_ >= 1 && maxPairDistance_ <= 0) return;

    auto distance = [&](const List<label>& group, bool physicalOnly) -> scalar
    {
        scalar maximum = 0;
        forAll(group, i) for (label j=i+1; j<group.size(); ++j)
        {
            const auto& a = eulerianFieldList[group[i]];
            const auto& b = eulerianFieldList[group[j]];
            scalar d = magSqr(a.position()-b.position());
            if (!physicalOnly)
            {
                d = physicalLocalization_ ? d/sqr(ri_) : 0;
                forAll(Xii_, k) d += sqr((a.XiR()[k]-b.XiR()[k])/Xii_[k]);
            }
            maximum = max(maximum, d);
        }
        return maximum;
    };
    std::stable_sort(pairs.begin(), pairs.end(), [&](const List<label>& a, const List<label>& b)
        { return distance(a, false) < distance(b, false); });
    const label budget = label(std::ceil(retainedPairFraction_*eulerianFieldList.size()));
    DynamicList<List<label>> accepted;
    label retained = 0;
    for (const auto& group : pairs)
    {
        if (maxPairDistance_ > 0 && distance(group, true) > sqr(maxPairDistance_)) continue;
        if (retained + group.size() > budget) continue;
        accepted.append(group);
        retained += group.size();
    }
    pairs = std::move(accepted);
}

template <class CloudType>
void Foam::mixParticleModel<CloudType>::KkdTreeLikeSearch
(
    const DynamicList<eulerianFieldData>& particleList,
    label l,
    label u,
    std::vector<label>& pInd,
    std::vector<label>& L,
    std::vector<label>& U    
) const
{
	//- Break the division if the particle list has length less than 2
    if (u - l <= 2)
    {
        //- Divide particles into groups of two or three
        L.push_back(l);

        U.push_back(u);

        return ;
    }

    label m = (l + u)/2;
    if ( (u - m) % 2 != 0 ) m++;

    auto iterL = pInd.begin();

    auto iterM = pInd.begin();

    auto iterU = pInd.begin();

    std::advance(iterL,l-1);

    std::advance(iterM,m-1);

    std::advance(iterU,u  );


    vector maxPosition(-GREAT, -GREAT, -GREAT);
    List<scalar> maxInXiR(Xii_.size(),-GREAT);

    vector minPosition(GREAT, GREAT, GREAT);
    List<scalar> minInXiR(Xii_.size(),GREAT);

    // Find minimum and maximum for each coordinate
    for (auto it = iterL; it != iterU; it++)
    {
        auto& pos = particleList[*it].position();
        for (label axis=0; axis<3; ++axis)
        {
            maxPosition[axis] = max(maxPosition[axis], pos[axis]);
            minPosition[axis] = min(minPosition[axis], pos[axis]);
        }

        forAll(Xii_,i)
        {
            maxInXiR[i] = std::max(maxInXiR[i],particleList[*it].XiR()[i]);
            minInXiR[i] = std::min(minInXiR[i],particleList[*it].XiR()[i]);
        }
    }

    //- Scaled/stretched distances between Max and Min in each direction
    //- Default is random mixing, overwritten if mixing distances greater than ri or fm
    scalar disMax = 0;
    label ncond = 0;

    if (physicalLocalization_)
        for (label axis=0; axis<3; ++axis)
        {
            const scalar span = (maxPosition[axis]-minPosition[axis])/ri_;
            if (span > disMax) { disMax = span; ncond = axis; }
        }
	scalar disXiR(0.0);
    forAll(Xii_,i)
    {
        disXiR = mag(maxInXiR[i] - minInXiR[i])/Xii_[i];
        if(disXiR > disMax)
        {
            disMax = disXiR;
            ncond = 3+i;
	    //Info << "ncond= " << ncond << endl;
            //ncond = particleList[*it].XiR()[i];
            //Info << "i: " << i << endl;
        }
    }

     // Tally which XiR axis won this split for the diagnostic histogram
    if (ncond >= 3 && (ncond - 3) < splitAxisHistogram_.size())
    {
        splitAxisHistogram_[ncond - 3]++;
    }

    lessArg comp(ncond);
    std::sort
    (
        iterL,
        iterU,
        [&](label& A, label& B) -> bool
        {
            if (comp(particleList[A],particleList[B])) return true;
            if (comp(particleList[B],particleList[A])) return false;
            return A < B;
        }
    );

    //- Recursive function calls for lower and upper branches of the particle list
    KkdTreeLikeSearch(particleList,l,m,pInd,L,U);

    KkdTreeLikeSearch(particleList,m+1,u,pInd,L,U);
};







template<class CloudType>
Foam::scalar Foam::mixParticleModel<CloudType>::readMixingConstant() const
{
    const dictionary& d = this->coeffDict();
    if (d.found("C_E"))
        WarningInFunction << "C_E is a deprecated alias for CE; use CE only." << nl;
    const scalar value = d.template lookupOrDefault<scalar>("CE",
        d.template lookupOrDefault<scalar>("C_E", 0.1));
    if (d.found("CE") && d.found("C_E") && value != readScalar(d.lookup("C_E")))
        FatalErrorInFunction << "Conflicting CE and C_E" << exit(FatalError);
    if (!(value > 0 && std::isfinite(value)))
        FatalErrorInFunction << "CE must be positive and finite" << exit(FatalError);
    return value;
}

template<class CloudType>
Foam::scalar Foam::mixParticleModel<CloudType>::pairMixingExtent
(
    const eulerianFieldData& p, const eulerianFieldData& q, scalar dt, scalar tau
) const
{
    if (!(tau > 0 && dt >= 0 && std::isfinite(tau) && std::isfinite(dt)))
        FatalErrorInFunction << "Invalid mixing time: dt=" << dt << ", tau=" << tau << exit(FatalError);
    minMixingTime_ = min(minMixingTime_, tau);
    maxMixingCourant_ = max(maxMixingCourant_, dt/tau);
    try
    {
        return secondConditioningNumerics::mixingExtent
            (dt, tau, mixingExtentModel_ == "modifiedCurl", p.Rand(), q.Rand());
    }
    catch (const std::exception& error)
    {
        FatalErrorInFunction << error.what() << "; dt=" << dt << ", tau=" << tau << exit(FatalError);
    }
    return 0;
}
