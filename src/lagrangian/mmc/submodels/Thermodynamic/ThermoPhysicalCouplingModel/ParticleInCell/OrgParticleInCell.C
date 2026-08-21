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

    forAll(mesh_.cells(),celli)
    {

        label superCelli = pManager.getSuperCellID(celli);

        DynamicList<particleType*> cellParticles = 
            pManager.getParticlesInSuperCell(this->owner(),superCelli);

        if (cellParticles.size() == 0)
        {
            // Do nothing. It will just use previous time-step value.
            // This if loop has to be included; otherwise hEqvETarget will be zero when Npc=0
        }
        else
        {
	//    if (this->owner().secondCondMixingEnabled() && cellParticles[iter]->secondCondFlag() != 1)
	//	continue;

            //- Find Average Temperature
            scalar totT = 0;
            scalar totWt = 0;

            forAll(cellParticles,iter)
            {
//		if (this->owner().secondCondMixingEnabled() && *cellParticles[iter]->secondCondFlag() != 1)
//	                continue;
                const particleType& p = *cellParticles[iter];       
                totWt += p.wt();
                totT += p.T() * p.wt(); 
            }

            TEqvETarget[celli] = totT/max(totWt,1e-15);


            forAll(YEqvETarget,specieI)
            {
//		if (this->owner().secondCondMixingEnabled() && cellParticles[iter]->secondCondFlag() != 1)
//			continue;

                if (!this->solveEqvSpecie()[specieI])
                {
                    YEqvETarget[specieI].primitiveFieldRef() = 0;
                } 
                else
                {
                    scalar totY = 0;
                    
                    forAll(cellParticles,iter)
                    {
                        const particleType& p = *cellParticles[iter];
                        totY += p.Y()[specieI] * p.wt();
                    }

                    YEqvETarget[specieI][celli] = totY/max(totWt,1e-15);
                }
            }
        }

    }
   
}

// ************************************************************************* //

