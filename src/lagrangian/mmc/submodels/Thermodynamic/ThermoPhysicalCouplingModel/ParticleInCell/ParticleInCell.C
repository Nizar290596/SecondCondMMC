/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParticleInCell.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleInCell<CloudType>::ParticleInCell
(
    const dictionary&,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    ThermoPhysicalCouplingModel<CloudType>(owner,Xi),
    mesh_(owner.mesh())
{
    this->Indicator() = 1.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleInCell<CloudType>::~ParticleInCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleInCell<CloudType>::EqvETargetValues
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{    
    // Reference to the particleNumber manager and super mesh
    const particleNumberController& pManager = this->owner().pManager();

    // When second-conditioning is active, restrict coupling to flagged
    // (subset) particles only. Cells without any flagged particle get
    // Indicator = 0 so the relaxation source is switched off there.
    const bool filterFlagged = this->owner().secondCondMixingEnabled();
    
    if (filterFlagged)
    {
        this->Indicator() = 0.0;
    }
    
    // Zero the target for species that are not solved (was done redundantly
    // inside the per-cell loop before; hoisted here so the whole field is
    // assigned once rather than N_cells times).
    forAll(YEqvETarget, specieI)
    {
        if (!this->solveEqvSpecie()[specieI])
        {
            YEqvETarget[specieI].primitiveFieldRef() = 0;
        }
    }
    
    // Single O(N_p) pass over the cloud: particles bucketed by super-cell.
    // Avoids the O(N_cells x N_p) scan that resulted from calling
    // getParticlesInSuperCell() once per cell.
    List<DynamicList<particleType*>> cellParticlesInSuperCell =
        pManager.getParticlesInSuperCellList(this->owner());

    const List<DynamicList<label>>& superCellCells =
        pManager.cellsInSuperCell();

    for
    (
        label superCelli = 0;
        superCelli < pManager.nSuperCells();
        ++superCelli
    )
    {
        const DynamicList<particleType*>& cellParticles =
            cellParticlesInSuperCell[superCelli];
            
            if (cellParticles.empty())
			{
				// Do nothing. It will just use previous time-step value.
				// This if loop has to be included; otherwise hEqvETarget will be zero when Npc=0
				// Keep previous-time-step target values in the member cells
				// (matches the original "Npc == 0" branch).
				continue;
			}
			
			// Single pass: accumulate T and solved species mean simultaneously.
			scalar totT  = 0;
			scalar totWt = 0;
			scalarField totY(YEqvETarget.size(), 0.0);

			for (particleType* pPtr : cellParticles)
			{
				if (filterFlagged && pPtr->secondCondFlag() != 1) continue;
				
				const scalar w = pPtr->wt();
				totWt += w;
				totT  += pPtr->T() * w;

				forAll(YEqvETarget, specieI)
				{
					if (this->solveEqvSpecie()[specieI])
					{
						totY[specieI] += pPtr->Y()[specieI] * w;
					}
				}
			}
			
			if (filterFlagged && totWt <= SMALL)
			{
				// No flagged particle in this super-cell - leave Indicator at 0
				continue;
			}
			
			const scalar invWt = 1.0/max(totWt, 1e-15);
			const scalar Tavg  = totT * invWt;
			
			const DynamicList<label>& memberCells = superCellCells[superCelli];
			for (const label celli : memberCells)
			{
				if (filterFlagged)
				{
					this->Indicator()[celli] = 1.0;
				}
				
				TEqvETarget[celli] = Tavg;

				forAll(YEqvETarget, specieI)
				{
					if (this->solveEqvSpecie()[specieI])
					{
						YEqvETarget[specieI][celli] = totY[specieI] * invWt;
					}
				}
			}
		}
}
// ************************************************************************* //

