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
    mixParticleModel<CloudType>(dict, owner, type, Xi, orderedRefVarNames(Xi))
{}


template <class CloudType>
Foam::mixParticleModel<CloudType>::mixParticleModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const mmcVarSet& Xi,
    const wordList& axisNames
)
:
    CloudMixingModel<CloudType>(dict,owner,type),

    XiR_(Xi),

    XiRNames_(orderedRefVarNames(Xi)),

    numXiR_(XiRNames_.size()),

    refAxisNames_(axisNames),

    Xii_(getXiNormalisation(axisNames)),

    ri_(readScalar(this->coeffDict().lookup("r_i"))),

    axisToRefVar_(mapAxesToRefVars(axisNames, XiRNames_)),

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
{}


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

    refAxisNames_(cm.refAxisNames_),

    Xii_(cm.Xii_),

    ri_(cm.ri_),

    axisToRefVar_(cm.axisToRefVar_),

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

    // Must match the primary constructor: the registered viscosity field is
    // the thermo one. A plain "mu" is created in createFields.H by copy and is
    // therefore NOT in the object registry, so looking it up here aborted any
    // clone() of a mixing model.
    mu_
    (
        this->owner().mesh().objectRegistry::lookupObject<volScalarField>
        ("thermo:mu")
    ),

    vb_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("vb")),
    

    pairingMethod_(cm.pairingMethod_),
    
    mixSubVolumes_(cm.mixSubVolumes_)

    //particlePairAlgorithm_(cm.particlePairAlgorithm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::mixParticleModel<CloudType>::buildParticleList()
{
    // ========================================================================
    // Build the particle list for this mixing stage.
    //
    // The three virtual hooks below are the only things a derived stage needs
    // to change. Everything else - the Eulerian interpolation, the parallel
    // exchange, the pairing - is shared, so a stage that pairs in its own
    // reference space still gets local, sub-volume and global pairing.
    //
    //   includeParticle()  which particles take part
    //   referenceVector()  the coordinates the k-d tree splits on
    //   reportPairing()    per-stage diagnostics, called at the end
    // ========================================================================

    // |grad Xi|^2 of every mmcVarSet reference variable, in XiRNames_ order.
    // Only used by the Cleary & Klimenko (aISO false) time scale.
    PtrList<volScalarField> magSqr_XiR_(XiRNames_.size());

    forAll(XiRNames_, refI)
    {
        magSqr_XiR_.set
        (
            refI,
            new volScalarField
            (
                magSqr(fvc::grad(this->XiR_.Vars(XiRNames_[refI]).field()))
            )
        );
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

    eulerianFieldDataList_.clear();

    particlePairs_.clear();

    const label nDims = this->nRefDims();

    scalarField xi(nDims, 0.0);

    // running index for particle position
    label particleInd=0;

    forAllIters(this->owner(), iter)
    {
        // Stage filter
        if (!this->includeParticle(iter()))
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

        // The eulerian data field also has to store the position for the
        // k-d tree later
        eulerianFields.position() = pos;

        // This stage's reference-space coordinates, written in refAxisNames_
        // order - the same order Xii_ was read in
        this->referenceVector(iter(), xi);

        eulerianFields.XiR() = xi;

        // Rand() is only consumed by the dense particle models, which build
        // their own lists; leave it defined but do not pay for an RNG here
        eulerianFields.Rand() = 0.0;

        eulerianFields.DEff() = DEff_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.D() = D_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.Dt() = Dt_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.DeltaE() = DeltaE_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.mu() = mu_intp_.interpolate(pos,cellI,faceI);

        eulerianFields.vb() = vb_intp_.interpolate(pos,cellI,faceI);

        // |grad Xi|^2 per pairing axis. Axes that are not mmcVarSet reference
        // variables have no Eulerian field and stay at zero.
        eulerianFields.magSqrRefVar().resize(nDims, 0.0);

        forAll(axisToRefVar_, axisI)
        {
            const label refI = axisToRefVar_[axisI];

            if (refI >= 0)
            {
                eulerianFields.magSqrRefVar()[axisI] =
                    magSqr_intp_[refI].interpolate(pos,cellI,faceI);
            }
        }

        eulerianFieldDataList_.append
        (
            std::move(eulerianFields)
        );

    }

    // Note: no early return here even when this process holds fewer than two
    // particles - correctParticleListParallel() is collective and every
    // process must reach it. findPairs() handles a short list itself.

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

    this->reportPairing();
}





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
        particleMixingProcessors_.size()*this->owner().size()
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
    // Assemble this stage's particle list and pair it
    buildParticleList();

    // Mix the list
    SmixList();
}


template <class CloudType>
void Foam::mixParticleModel<CloudType>::referenceVector
(
    const particleType& p,
    scalarField& xi
) const
{
    // Default pairing space: the mmcVarSet reference variables. XiRNames_ is
    // ordered by rVarInXiR() index, so component i of xi and Xii_[i] both
    // refer to the same reference variable.
    forAll(xi, i)
    {
        xi[i] = p.XiR()[i];
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

    // Sub-volume construction is only defined for a one-dimensional pairing
    // space (mixingSubVolumes takes a single normalisation factor). Test the
    // dimension of THIS stage's space, not the mmcVarSet reference count.
    if (Xii_.size() > 1)
        FatalError << "Mixing model " << this->modelType() << " pairs in a "
            << Xii_.size() << "-dimensional reference space ("
            << refAxisNames_ << ")." << nl
            << "pairingMethod subVolumes requires a single axis; use "
            << "local or global instead."
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
Foam::wordList
Foam::mixParticleModel<CloudType>::orderedRefVarNames(const mmcVarSet& Xi)
{
    // rVarInXiR() maps a reference variable name to the index it occupies in
    // the particle's XiR() array. Invert it, so that names[i] names XiR()[i].
    //
    // Do NOT use rVarInXi().toc() for this: HashTable::toc() returns keys in
    // hash order, which has nothing to do with the index rVarInXiR() assigned,
    // so any list built that way can pair a normalisation factor - or a
    // |grad Xi|^2 field - with the wrong reference variable.
    const HashTable<label, word>& refIndex = Xi.rVarInXiR();

    wordList names(refIndex.size());

    forAllConstIters(refIndex, iter)
    {
        const label i = iter.val();

        if (i < 0 || i >= names.size())
        {
            FatalErrorInFunction
                << "Reference variable " << iter.key() << " has index " << i
                << " which is out of range for " << names.size()
                << " reference variables"
                << exit(FatalError);
        }

        names[i] = iter.key();
    }

    return names;
}


template <class CloudType>
Foam::labelList Foam::mixParticleModel<CloudType>::mapAxesToRefVars
(
    const wordList& axisNames,
    const wordList& refVarNames
)
{
    labelList map(axisNames.size(), -1);

    forAll(axisNames, axisI)
    {
        forAll(refVarNames, refI)
        {
            if (axisNames[axisI] == refVarNames[refI])
            {
                map[axisI] = refI;
                break;
            }
        }
    }

    return map;
}


template <class CloudType>
List<scalar> Foam::mixParticleModel<CloudType>::getXiNormalisation
(
    const wordList& axisNames
) const
{
    // Normalisation factors for this stage's own pairing axes, read from this
    // stage's own coeffDict - so each mixing stage is configured entirely
    // within its own <model>Coeffs sub-dictionary.
    const dictionary XiDict(this->coeffDict().subDict("Xim_i"));

    Info << nl << "The Xim_i parameters of " << this->modelType()
         << " are: " << XiDict << endl;

    if (axisNames.empty())
    {
        FatalErrorInFunction
            << "Mixing model " << this->modelType() << " has no pairing axes."
            << nl
            << "At least one reference variable, or an explicit referenceAxes "
            << "entry, is required."
            << exit(FatalError);
    }

    List<scalar> Xii(axisNames.size());

    // Indexed by position in axisNames, which is also the order that
    // referenceVector() writes its components - so Xii[i] always normalises
    // component i and the two can never drift apart.
    forAll(axisNames, i)
    {
        const word key(axisNames[i] + "_m");

        if (!XiDict.found(key))
        {
            FatalErrorInFunction
                << "Missing normalisation factor " << key << " in "
                << this->modelType() << "Coeffs/Xim_i" << nl
                << "One entry is required per pairing axis. This model pairs "
                << "on: " << axisNames
                << exit(FatalError);
        }

        Xii[i] = readScalar(XiDict.lookup(key));
    }

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

    // A list of fewer than two particles has no pair to form. Without this the
    // loop below still reads pInd[0] and pInd[1]: KkdTreeLikeSearch() takes its
    // base-case branch immediately and pushes L=1, U=size, which the consumer
    // reads as the pair (pInd[0], pInd[1]) regardless of how long pInd is.
    // Empty and single-particle lists are routine for a stage that pairs a
    // subset of the cloud on one process.
    //
    // This guard belongs here rather than in buildParticleList(): returning
    // early from there would skip correctParticleListParallel(), which is
    // collective, and deadlock every other process.
    if (eulerianFieldList.size() < 2)
        return;

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


    //scalar maxInX = -GREAT;
    //scalar maxInY = -GREAT;
    //scalar maxInZ = -GREAT;
    List<scalar> maxInXiR(Xii_.size(),-GREAT);

    //scalar minInX = GREAT;
    //scalar minInY = GREAT;
    //scalar minInZ = GREAT;
    List<scalar> minInXiR(Xii_.size(),GREAT);

    // Find minimum and maximum for each coordinate
    for (auto it = iterL; it != iterU; it++)
    {
        auto& pos = particleList[*it].position();
        //maxInX = std::max(maxInX,pos.x());
        //maxInY = std::max(maxInY,pos.y());
        //maxInZ = std::max(maxInZ,pos.z());

        //minInX = std::min(minInX,pos.x());
        //minInY = std::min(minInY,pos.y());
        //minInZ = std::min(minInZ,pos.z());

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

    //scalar disX = (maxInX - minInX)/ri_;
    //if(disX > disMax)
    //{
        //disMax = disX;
        //ncond = 0;
    //}

    //scalar disY = (maxInY - minInY)/ri_;
    //if(disY > disMax)
    //{
        //disMax = disY;
        //ncond = 1;
    //}

    //scalar disZ = (maxInZ - minInZ)/ri_;
    //if(disZ > disMax)
    //{
        //disMax = disZ;
        //ncond = 2;
    //}
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
            return comp(particleList[A],particleList[B]);
        }
    );

    //- Recursive function calls for lower and upper branches of the particle list
    KkdTreeLikeSearch(particleList,l,m,pInd,L,U);

    KkdTreeLikeSearch(particleList,m+1,u,pInd,L,U);
};





