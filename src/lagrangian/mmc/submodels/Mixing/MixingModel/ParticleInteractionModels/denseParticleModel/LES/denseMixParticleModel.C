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

#include "denseMixParticleModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::denseMixParticleModel<CloudType>::denseMixParticleModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudMixingModel<CloudType>(dict,owner,type),
  
    DEff_(owner.mesh().objectRegistry::lookupObject<volScalarField>("DEff")),

    particleList_(owner.mesh().C().size()),

    eulerianFieldsList_(owner.mesh().C().size())
{
    Info << "Note: you are using a dense particle mixing model!" << endl;
    Info << "The number of supercells should equal the number of LES cells." << endl;
    Info << "Uses (non-interpolated) D_eff of LES cells. Only correct if n_LES = n_superCell" << endl;
    
    // Check that the mesh and super mesh match
    if 
    (
        this->owner().pManager().nSuperCellsOnProcessor() 
     != this->owner().mesh().C().size()
    )
    FatalErrorInFunction 
        << "SuperMesh and LES mesh do not have the same cell number"
        << exit(FatalError);
    
}


template <class CloudType>
Foam::denseMixParticleModel<CloudType>::denseMixParticleModel
(
    const denseMixParticleModel<CloudType>& cm
)
:
    CloudMixingModel<CloudType>(cm),
    
    DEff_(this->owner().mesh().objectRegistry::lookupObject<volScalarField>("DEff")),
    
    particleList_(this->owner().mesh().C().size()),

    eulerianFieldsList_(this->owner().mesh().C().size())
{
    Info << "Note: you are using a dense particle mixing model!" << endl;
    Info << "The number of supercells should equal the number of LES cells." << endl;
    Info << "Uses (non-interpolated) D_eff of LES cells. Only correct if n_LES = n_superCell" << endl;
    
    // Check that the mesh and super mesh match
    if 
    (
        this->owner().pManager().nSuperCellsOnProcessor() 
     != this->owner().mesh().C().size()
    )
    FatalErrorInFunction 
        << "SuperMesh and LES mesh do not have the same cell number"
        << exit(FatalError);    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::denseMixParticleModel<CloudType>::~denseMixParticleModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::denseMixParticleModel<CloudType>::buildParticleList()
{
    interpolationCellPoint<scalar> DEff_intp_(this->DEff_);

    StochasticLib1 rand(time(0));
    
    // =========================================================================
    // Note: Super mesh and LES mesh match -- iterartes only over the LES mesh
    // =========================================================================

    // Clear old lists
    forAll(particleList_,celli)
    {
        particleList_[celli].clear();
        eulerianFieldsList_[celli].clear();
    }


    // create mixParticle and append to list
    forAllIters(this->owner(), iter)
    {        
        eulerianFieldParticleValues eulFields;
        
        eulFields.Rand() = rand.Random();
                
        // Because super mesh and LES mesh match iter().cell() == superCellI
        const label celli = iter().cell();
        
        eulFields.sqrFilter() = pow(this->owner().mesh().V()[celli],2./3.);

        eulFields.DEff() = DEff_[celli]; 
 
        eulerianFieldsList_[celli].append(std::move(eulFields));
        particleList_[celli].append(iter.get());
    }
} 

template<class CloudType>
void Foam::denseMixParticleModel<CloudType>::Smix()
{
    buildParticleList();
    
    // Loop over all cells and mix particles in each cell with each other
    forAll(particleList_, celli)
    {
        SmixList(particleList_[celli],eulerianFieldsList_[celli]);
    }
}

// ************************************************************************* //
