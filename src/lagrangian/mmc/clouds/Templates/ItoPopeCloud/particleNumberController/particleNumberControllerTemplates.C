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

\*---------------------------------------------------------------------------*/

template<class CloudType>
void Foam::particleNumberController::correct
(
    CloudType& cloud
)
{
    // If number control is switched off - do nothing
    if (!numCtrlOn_)
        return;
    
    using ParticleType = typename CloudType::particleType;
   
    // get particles of neighboring processors
    auto particlesInSuperCell = collectParticleList(cloud);
   
    List<int> seeds = synchronizeRandomSeed();
    
    forAll(particlesInSuperCell,superCellI)
    {
        StochasticLib1 rndGen(seeds[superCellI]);
        
        // Check if this super cell is managed by this processor
        if (superCellOnProcessor_[superCellI])
        {   
            auto& cellParticles(particlesInSuperCell[superCellI]);
            
            const label nc = cellParticles.size();

            //- Particles have their weight doubled and then are deleted with
            //- a probability of 50% 
            
            if (nc > Nhi_)
            {
                std::sort(cellParticles.begin(),cellParticles.end(),lessWt<ParticleType>());
            
                const label particlesToDelete = nc-Npc_;
                label particlesDeleted = 0;
                label i = 0;
                while (particlesDeleted < particlesToDelete && i < nc)
                {
                    // Double the weight of the particle
                    if (cellParticles[i].second() != nullptr)
                    {
                        cellParticles[i].second()->wt() = 
                            2.0*cellParticles[i].second()->wt();
                        cellParticles[i].second()->m() = 
                            2.0*cellParticles[i].second()->m();
                    }

                    // Set the flag if the particle would be deleted
                    bool deleteParticle = false;
                    if (rndGen.Random() >= 0.5)
                    {
                        deleteParticle = true;
                        particlesDeleted++;
                    }

                    // Can only perform clone and delete on local particle set
                    if (deleteParticle && cellParticles[i].second() != nullptr)
                        cloud.deleteParticle(*(cellParticles[i].second()));

                    i++;   
                }
            }
            else if (nc != 0 && nc < Nlo_)
            {
                std::sort(cellParticles.begin(),cellParticles.end(),lessWt<ParticleType>());
            
                label nclone = Npc_ - nc;
            
                nclone = min(nclone, nc);
            
                for(label i = nc-1; i > nc-nclone-1; i--)
                {
                    if (cellParticles[i].second() != nullptr)
                    {
                        cellParticles[i].second()->wt() = 
                            0.5*cellParticles[i].second()->wt();
                        cellParticles[i].second()->m() = 
                            0.5*cellParticles[i].second()->m();
                
                        const ParticleType* p = cellParticles[i].second();
                    
                        cloud.append(p->clone().ptr());
                    }
                }
            }
        }
    }    
}


template<class CloudType>
Foam::List
<
    Foam::DynamicList
    <
        Foam::Tuple2
        <
            Foam::scalar,
            typename CloudType::particleType*    
        >
    > 
>
Foam::particleNumberController::collectParticleList(CloudType& cloud)
{
    using ParticleType = typename CloudType::particleType;    
    
    List<DynamicList<scalar>> particleInSuperCellWeights(nSuperCells_);
    
    // Store the particle weights, required for the sorting later
    for (auto it=cloud.begin(); it != cloud.end(); ++it)
    {
        it().sCell() = superCellForCell_[it().cell()];
        const label superCellI = superCellForCell_[it().cell()];

        particleInSuperCellWeights[superCellI].append(it().wt());
    }

    if (Pstream::parRun())
    {
        // Send this to all procesors that also contain this super cell
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // To which processors does the list need to be send
        forAll(particleInSuperCellWeights, superCellI)
        {
            if (superCellOnProcessor_[superCellI])
            {
                // To which processors does the list need to be send
                for (label procI : superCellToProcID_[superCellI])
                {
                    if (procI != Pstream::myProcNo())
                    {
                        UOPstream toBuffer(procI,pBufs);
                        toBuffer << particleInSuperCellWeights[superCellI];
                    }
                }
            }
        }

        pBufs.finishedSends();
        
        forAll(particleInSuperCellWeights, superCellI)
        {
            if (superCellOnProcessor_[superCellI])
            {
                // Receive from processor procI
                for (label procI : superCellToProcID_[superCellI])
                {
                    if (procI != Pstream::myProcNo())
                    {
                        List<scalar> particleWeightsOfProcI;
                        
                        UIPstream fromBuffer(procI,pBufs);
                        fromBuffer >> particleWeightsOfProcI;
                        particleInSuperCellWeights[superCellI].append
                        (
                            std::move(particleWeightsOfProcI)
                        );
                    }
                }
            }
        }
    }
    

    List
    <
        DynamicList
        <
            Tuple2
            <
                scalar,
                ParticleType*    
            >
        >
    >  particlesInSuperCell(nSuperCells_);
    
    List<label> numberOfLocalParticlesInSuperCell(nSuperCells_,0);
    
    for (auto it=cloud.begin(); it != cloud.end(); ++it)
    {
        const label superCellI = superCellForCell_[it().cell()];
        particlesInSuperCell[superCellI].append
        (
            Tuple2<scalar,ParticleType*>(it().wt(),it.get())
        );
        numberOfLocalParticlesInSuperCell[superCellI]++;
    }
    
    // Now add the processor particles with nullptr 
    forAll(particleInSuperCellWeights,superCellI)
    {
        forAll(particleInSuperCellWeights[superCellI],i)
        {
            if (i >= numberOfLocalParticlesInSuperCell[superCellI])
            {
                particlesInSuperCell[superCellI].append
                (
                    Tuple2<scalar,ParticleType*>
                    (
                        particleInSuperCellWeights[superCellI][i],
                        nullptr
                    )
                );                        
            }
        }
    }
    
    return particlesInSuperCell;
}


template<class FieldType>
Foam::List<Foam::DynamicList<typename FieldType::value_type>>
Foam::particleNumberController::fieldDataInSuperCell
(
    const FieldType& field
) const
{   
    typedef typename FieldType::value_type Type;
    
    // Local list used for sending
    List<DynamicList<Type>> localValuesInSuperCell(nSuperCells_);
    
    forAll(localValuesInSuperCell,superCellI)
    {
        for (label cellI : cellsInSuperCell_[superCellI])
            localValuesInSuperCell[superCellI].append(field[cellI]);
    }
    
    // List to return
    List<DynamicList<Type>> valuesInSuperCell(nSuperCells_);
    
    if (Pstream::parRun())
    {
        // Send this local cell positions
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        
        // Note to order of sending information with Pstream:
        // As all processors read the same super mesh the loop over all 
        // super cells is the same on each processor, hence the order of 
        // writing and reading into the Pstream buffer is the same.
        forAll(localValuesInSuperCell, superCellI)
        {
            if (superCellOnProcessor_[superCellI])
            {
                // To which processors does the list need to be send
                for (label procI : superCellToProcID_[superCellI])
                {
                    if (procI != Pstream::myProcNo())
                    {
                        UOPstream toBuffer(procI,pBufs);
                        toBuffer << localValuesInSuperCell[superCellI];
                    }
                }
            }
        }
        
        pBufs.finishedSends();
        
        forAll(valuesInSuperCell, superCellI)
        {
            if (superCellOnProcessor_[superCellI])
            {
                for (label procI : superCellToProcID_[superCellI])
                {
                    // Keep the order for all processors the same
                    if (procI != Pstream::myProcNo())
                    {
                        DynamicList<Type> temp;
                        UIPstream fromBuffer(procI,pBufs);
                        fromBuffer >> temp;
                        valuesInSuperCell[superCellI].append
                        (
                            std::move(temp)
                        );
                    }
                    else
                    {
                        valuesInSuperCell[superCellI].append
                        (
                            std::move(localValuesInSuperCell[superCellI])
                        );
                    }
                }
            }
        }
    }
    else
    {
        // If not parallel 
        valuesInSuperCell = std::move(localValuesInSuperCell);
    }
    
    return valuesInSuperCell;
}


template<class CloudType>
Foam::List<Foam::DynamicList<typename CloudType::particleType*>>
Foam::particleNumberController::getParticlesInSuperCellList
(
    CloudType& cloud
) const
{
    typedef typename CloudType::particleType ParticleType;
    List<DynamicList<ParticleType*>> pList(nSuperCells_);
    
    // Loop over all particles and check in which super cell they are 
    for (auto it=cloud.begin(); it != cloud.end(); ++it)
    {
        pList[getSuperCellID(it().cell())].append(it.get());
    }
    
    return pList;
}


template<class CloudType>
Foam::DynamicList<typename CloudType::particleType*>
Foam::particleNumberController::getParticlesInSuperCell
(
    CloudType& cloud,
    const label superCellI
) const
{
    typedef typename CloudType::particleType ParticleType;
    DynamicList<ParticleType*> pList;
    
    // Loop over all particles and check in which super cell they are 
    for (auto it=cloud.begin(); it != cloud.end(); ++it)
    {
        // if particle is in super cell superCellI add to list
        if (getSuperCellID(it().cell()) == superCellI)
            pList.append(it.get());
    }
    
    return pList;
}
