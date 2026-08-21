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

\*---------------------------------------------------------------------------*/

#include "premixedMixParticleModel.H"
#include <unordered_set>
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::premixedMixParticleModel<CloudType>::premixedMixParticleModel
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
    
    ri_(readScalar(this->coeffDict().lookup("ri"))),

    Xii_(getXiNormalisation()),

    D_(owner.mesh().objectRegistry::lookupObject<volScalarField>("D")),

    Dt_(owner.mesh().objectRegistry::lookupObject<volScalarField>("Dt")),

    DeltaE_(owner.mesh().objectRegistry::lookupObject<volScalarField>("DeltaE")),

    mu_(owner.mesh().objectRegistry::lookupObject<volScalarField>("thermo:mu")),
    
    TF_(owner.mesh().objectRegistry::lookupObject<volScalarField>("TF")),
    
    ATFkd_(this->coeffDict().lookup("ATFkd")),
    
    pairingMethod_(this->coeffDict()),

    mixSubVolumes_(owner.mesh(),this->coeffDict(),ri_,Xii_[0])
{}


template <class CloudType>
Foam::premixedMixParticleModel<CloudType>::premixedMixParticleModel
(
    const premixedMixParticleModel<CloudType>& cm
)
:
    CloudMixingModel<CloudType>(cm),
    
    XiR_(cm.XiR_),

    XiRNames_(cm.XiRNames_),
    
    numXiR_(XiRNames_.size()),

    ri_(readScalar(this->coeffDict().lookup("ri"))),

    Xii_(getXiNormalisation()), 

    //initialise the fields needed for aISO
    D_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("D")),

    Dt_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("Dt")),

    DeltaE_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("DeltaE")),

    mu_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("mu")),
    
    TF_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("TF")),
    
    ATFkd_(this->coeffDict().lookup("ATFkd")),
    
    pairingMethod_(cm.pairingMethod_),
    
    mixSubVolumes_(cm.mixSubVolumes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::premixedMixParticleModel<CloudType>::buildParticleList()
{
    // ========================================================================
    // Build local particleList

    PtrList<volScalarField> magSqr_XiR_(XiRNames_.size());
    
    label II = 0;
    
    for (const word& nameI : XiRNames_)
    {
        magSqr_XiR_.set
        (
            II,
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

    //for aISO
    interpolationCellPoint<scalar> D_intp_(this->D_);

    interpolationCellPoint<scalar> Dt_intp_(this->Dt_);

    interpolationCellPoint<scalar> DeltaE_intp_(this->DeltaE_);

    interpolationCellPoint<scalar> mu_intp_(this->mu_);
    
    interpolationCellPoint<scalar> TF_intp_(this->TF_);
  
    // clear particle list from old data
    particleList_.clear();
    
    StochasticLib1 rand(time(0));
    
    eulerianFieldDataList_.clear();
    
    // running index for particle position
    label particleInd=0;
    
    forAllIters(this->owner(), iter)
    {
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
        
        eulerianFields.D() = D_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.Dt() = Dt_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.DEff() = eulerianFields.D() + eulerianFields.Dt();
            
        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos,cellI,faceI);
            
        eulerianFields.mu() = mu_intp_.interpolate(pos,cellI,faceI);
        
        eulerianFields.TF() = TF_intp_.interpolate(pos,cellI,faceI);
               
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
    
    // if run in parallel get all required particles of neighbouring processors
    if 
    (
            Pstream::parRun() 
         && pairingMethod_.method() != particlePairingMethod::localPairing
    )
    {
        // findPairs is called in correctParticleListParallel
        correctParticleListParallel();
    }
    else
    {
        findPairs(eulerianFieldDataList_,particlePairs_);
    }
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::correctParticleListParallel()
{  
    particleMixingProcessors_ = getParticleMixingProcessors();

    // Collect the eulerian data from other processors
    collectEulerianDataFields();

    // Find the particle pairs to mix
    findPairs(eulerianFieldDataList_,particlePairs_);
  
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
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::collectEulerianDataFields()
{
    // First send all the eulerian data fields and create the pairing lists
    // Then only send pairs that mix with particles located on the current
    // processor 
    
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // send the list of particles of this processor to neighbour 
    // processors
    for (auto& procI : particleMixingProcessors_)
    {
        if (procI != Pstream::myProcNo())
        {
            UOPstream toBuffer(procI,pBufs);
            toBuffer << eulerianFieldDataList_.size();
            for (auto& e : eulerianFieldDataList_)
                toBuffer << e;
        }
    }
    
    pBufs.finishedSends();

    startIndexOfParticle_.clear();
    startIndexOfParticle_.resize(Pstream::nProcs(),-1);
    label previousParticleSize = 0;

    // Loop over all processors to update the startIndexOfParticle list
    // it is important that all processors have the same order of Eulerian
    // fields. Otherwise they might calculate different pairings!

    DynamicList<eulerianFieldData> tlocalEulerianFields = std::move(eulerianFieldDataList_);
    eulerianFieldDataList_.clear();

    // Estimate space for eulerianFields
    eulerianFieldDataList_.reserve
    (
        particleMixingProcessors_.size()*this->owner().size()
    );

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
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::collectParticleData()
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
    const label numLocalParticles = this->owner().size();

    // Reserve some space
    forAll(particlesToSendToProcessor,i)
    {
        particlesToSendToProcessor[i].reserve(0.05*numLocalParticles);
    }
    
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

                eulerianFieldDataList_[startIndexOfParticle_[procI]+procParticleInd].particleIndex() = particleInd++;
            }
        }
    }
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::Smix()
{
    // Build the particle list
    buildParticleList();
    
    // Mix the list
    SmixList();
}


template<class CloudType>
Foam::List<Foam::label>
Foam::premixedMixParticleModel<CloudType>::getParticleMixingProcessors()
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

    // Only works for one reference variable 
    if (numXiR_ > 1)
        FatalError << "More than one reference variable selected."<<nl
            << "Only particle pairing local and global are possible"
            << exit(FatalError);

    
    // get the sub-volume
    return mixSubVolumes_.getSubVolume(Pstream::myProcNo());
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::SmixList()
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

            // Only call mixpair for processor local particles, as remote
            // particles are mixed on their respective processor

            if (e1.local() || e2.local())
                mixpair
                (
                    particleList_[e1.particleIndex()],e1,
                    particleList_[e2.particleIndex()],e2,
                    deltaT
                );
            if (e2.local() || e3.local())
                mixpair
                (
                    particleList_[e2.particleIndex()],e2,
                    particleList_[e3.particleIndex()],e3,
                    deltaT
                );
        }
    }
}


template <class CloudType>
List<scalar> Foam::premixedMixParticleModel<CloudType>::getXiNormalisation() 
{
    // dictionary to read the normalisation parameters for 
    // the reference variables 
    const dictionary XiDict(this->coeffDict().subDict("Ximi"));

    Info << nl << "The Ximi parameters are: "<< XiDict << endl;

    List<scalar> Xii(numXiR_);

    label i=0;
    for (const word& refVarName :this->XiRNames())
    {
        Xii[i++] = readScalar(XiDict.lookup(refVarName+"m"));
    }

    return Xii;
}


template<class CloudType>
void Foam::premixedMixParticleModel<CloudType>::findPairs
(
    const DynamicList<eulerianFieldData>& eulerianFieldList,
    DynamicList<List<label>>& pairs
) const
{
    // Clear particle pairs first
    pairs.clear();

    // Keeping track of indices for premixedkdTreeLikeSearch
    std::vector<label> L;
    std::vector<label> U;
    L.reserve(eulerianFieldList.size());
    U.reserve(eulerianFieldList.size());
    
    // create an index list for the particle data
    std::vector<label> pInd(eulerianFieldList.size());
    std::iota(pInd.begin(),pInd.end(),0);
    
    premixedkdTreeLikeSearch(eulerianFieldList,1,eulerianFieldList.size(),pInd,L,U);

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
    }
}

// ************************************************************************* //

template <class CloudType>
void Foam::premixedMixParticleModel<CloudType>::premixedkdTreeLikeSearch
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
    
    // Only works for one reference variable 
    if (numXiR_ > 1)
        FatalError << "More than one reference variable selected."<<nl
                   << "Premixed kdTree search is implemented only for a single reference variable."
                   << exit(FatalError);

    scalar maxInX = -GREAT;
    scalar maxInY = -GREAT;
    scalar maxInZ = -GREAT;
    scalar maxInXiR = -GREAT;
    
    scalar minInX = GREAT;
    scalar minInY = GREAT;
    scalar minInZ = GREAT;
    scalar minInXiR = GREAT;
    
    // Thickening factors of the particles with max and min reference variable values
    scalar maxInXiRTF = -GREAT;
    scalar minInXiRTF = GREAT;
    
    // Maximum thickening factor of all particles in the list
    scalar maxTF = -GREAT;
    
    // Find minimum and maximum for each coordinate
    for (auto it = iterL; it != iterU; it++)
    {
        auto& pos = particleList[*it].position();
        maxInX = std::max(maxInX,pos.x());
        maxInY = std::max(maxInY,pos.y());
        maxInZ = std::max(maxInZ,pos.z());
        
        minInX = std::min(minInX,pos.x());
        minInY = std::min(minInY,pos.y());
        minInZ = std::min(minInZ,pos.z());
        
        maxTF = std::max(maxTF,particleList[*it].TF());
        
        // Attention: This part only works for one reference variable!             
        if (particleList[*it].XiR()[0] > maxInXiR)
        {
            maxInXiR = particleList[*it].XiR()[0];
            maxInXiRTF = particleList[*it].TF();
        }    
        if (particleList[*it].XiR()[0] < minInXiR)
        {
            minInXiR = particleList[*it].XiR()[0];
            minInXiRTF = particleList[*it].TF();
        }  
    }

    //- Scaled/stretched distances between Max and Min in each direction
    //- Default is random mixing, overwritten if mixing distances greater than ri or fm
    scalar disMax = 0;
    label ncond = 0;

    scalar disX = (maxInX - minInX)/ri_;
    if(disX > disMax)
    {
        disMax = disX;
        ncond = 0;
    }

    scalar disY = (maxInY - minInY)/ri_;
    if(disY > disMax)
    {
        disMax = disY;
        ncond = 1;
    }

    scalar disZ = (maxInZ - minInZ)/ri_;
    if(disZ > disMax)
    {
        disMax = disZ;
        ncond = 2;
    }

    // - When ATF is used the reference variable gradients are reduced.
    //   Therefore, characteristic distance in the reference space should proportionally decrease thus preserving localness of mixing in composition.
    scalar disXiR = 0.0;
    if (ATFkd_)
    {
        if((maxInXiR - 0.5) * (minInXiR - 0.5) <= 0.0) // If particles are on the different sides from the flame center (XiR = 0.5) use overall max thickening factor 
        {
            disXiR = (maxInXiR - minInXiR) / (Xii_[0]/maxTF);
        }
        else // Use max thickening factor chosen between the particles with max and min reference variable values
        {
            scalar XiRTF = std::max(maxInXiRTF, minInXiRTF);
            disXiR = (maxInXiR - minInXiR) / (Xii_[0]/XiRTF);
        }
    }
    else 
    {
        disXiR = (maxInXiR - minInXiR)/Xii_[0];    
    }
    if(disXiR > disMax)
    {
        disMax = disXiR;
        ncond = 3;
    }


    lessArg comp(ncond);
    std::sort
    (
        iterL,
        iterU,
        [&](label& A, label& B) -> bool
        {
            return comp(particleList[A],particleList[B]);
        }
    );

    //- Recursive function calls for lower and upper branches of the particle list
    premixedkdTreeLikeSearch(particleList,l,m,pInd,L,U);

    premixedkdTreeLikeSearch(particleList,m+1,u,pInd,L,U);
};




