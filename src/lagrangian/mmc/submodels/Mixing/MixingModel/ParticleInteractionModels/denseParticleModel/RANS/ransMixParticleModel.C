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

#include "ransMixParticleModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ransMixParticleModel<CloudType>::ransMixParticleModel
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
  
    tauTurb_(owner.mesh().objectRegistry::lookupObject<volScalarField>("tauTurb")),

    particleList_(owner.pManager().superMesh().nCells()),
    
    eulerianFieldsList_(owner.pManager().superMesh().nCells())    
{}


template <class CloudType>
Foam::ransMixParticleModel<CloudType>::ransMixParticleModel
(
    const ransMixParticleModel<CloudType>& cm
)
:
    CloudMixingModel<CloudType>(cm),

    XiR_(cm.XiR_),

    XiRNames_(cm.XiRNames_),
    
    numXiR_(XiRNames_.size()),
    
    tauTurb_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("tauTurb")),

    particleList_(this->owner().pManager().superMesh().nCells()),
    
    eulerianFieldsList_(this->owner().pManager().superMesh().nCells())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ransMixParticleModel<CloudType>::~ransMixParticleModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::ransMixParticleModel<CloudType>::buildParticleList()
{
    // Reference to particle manager
    const particleNumberController& pManager = this->owner().pManager();
  
    //- clear particle list
    forAll(particleList_, i)
    {
        particleList_[i].clear();
        eulerianFieldsList_[i].clear();
    }


    StochasticLib1 rand(time(0));

    // get value tauTurb in each super cell
    List<DynamicList<scalar>> tauTurbInSuperCell = 
        pManager.fieldDataInSuperCell(tauTurb_);


    // Calculate average turbulent time scale for a super cell    
    List<scalar> tauTurbSuperCell_(pManager.superMesh().nCells(), 0.);
    List<label> nCellsInSuperCell_(pManager.superMesh().nCells(), 0);

    // Loop over all super cells
    forAll(tauTurbInSuperCell,superCelli)
    {
        // Loop over all cells in the super cell
        forAll(tauTurbInSuperCell[superCelli],i)
        {
            tauTurbSuperCell_[superCelli] += tauTurbInSuperCell[superCelli][i];
            nCellsInSuperCell_[superCelli]++;
        }
        tauTurbSuperCell_[superCelli] /= (nCellsInSuperCell_[superCelli] + VSMALL);
    }

    //- create mixParticle and append to list
    forAllIters(this->owner(), iter)
    {
        eulerianFieldData eulField; 
        
        // The eulerian data field also has to store the position for the 
        // k-d tree later
        eulField.position() = iter().position();
        
        // Also store the reference variable
        eulField.XiR() = iter().XiR();
        
        eulField.Rand() = rand.Random();
                
        label superCelli = pManager.getSuperCellID(iter().cell());
        
        eulField.tauTurb() = tauTurbSuperCell_[superCelli];

        particleList_[superCelli].append(iter.get());
        
        // Store index of particle
        eulField.particleIndex() = particleList_[superCelli].size()-1;
        eulerianFieldsList_[superCelli].append(std::move(eulField));
    }
} 


template<class CloudType>
void Foam::ransMixParticleModel<CloudType>::Smix()
{
    buildParticleList();
    
    //- mix in every supercell
    forAll(particleList_, i)
    {
        SmixList(particleList_[i],eulerianFieldsList_[i]);   
    }
}


// ************************************************************************* //
