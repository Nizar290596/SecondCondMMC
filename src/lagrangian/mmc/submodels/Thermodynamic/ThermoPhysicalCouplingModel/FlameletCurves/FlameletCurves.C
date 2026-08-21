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

#include "FlameletCurves.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FlameletCurves<CloudType>::FlameletCurves
(
    const dictionary&,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    ThermoPhysicalCouplingModel<CloudType>(owner,Xi),
    primarySpecies_("undefined"),
    independentYEqvCurve_(false),
    TEqvTargetCurve_(owner.pManager().nSuperCells()),
    YEqvTargetCurve_(owner.pManager().nSuperCells())
{
    // Read in which species are to be used in the coupling, and the primary
    // species for dynamic matching
    const dictionary couplingDict
    (
        this->owner().cloudProperties().subDict("thermophysicalCoupling")
    );
    
    
    // Primary species is the species used for selecting the equivalent species
    // target curve
    primarySpecies_ = couplingDict.get<word>("primarySpecies");
    
    // If primary species is temperature then species target curves
    // are not independent of the temperature target curve
    if (primarySpecies_ == "T")
        independentYEqvCurve_ = false;
    else
        independentYEqvCurve_ = true;
    
    this->Indicator() = 1.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FlameletCurves<CloudType>::~FlameletCurves()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FlameletCurves<CloudType>::findEquivalentTempSpeciesCurves()
{
    typedef typename CloudType::particleType particleType;
    
    
    // Reference to particle management
    auto& pManager = this->owner().pManager();
    
    // Create interpolation tables
    interpolationLookUpTable<scalar> TEqvInterpTab
    (
        "TEqv",
        this->owner().mesh().time().constant(),
        this->owner().mesh()
    );

    // Get the list of particles in each superCell
    List<DynamicList<particleType*>> cellParticlesInSuperCell = 
        pManager.getParticlesInSuperCellList(this->owner());
        
    // Find curves in each supercell
    for (label superCelli=0; superCelli < pManager.nSuperCells(); superCelli++)
    {        
        DynamicList<particleType*> cellParticles 
            = cellParticlesInSuperCell[superCelli];

        //- Find equivalent temperature curve
        scalar TEqvTmp;

        for (label i = 1; i < TEqvInterpTab.output().size(); i++)
        {
            scalar totT = 0;

            for (particleType* pPtr : cellParticles)
                totT += sqr
                    (
                        pPtr->T() 
                      - TEqvInterpTab.lookUp(pPtr->XiC(this->cVarName()))[i]
                    ) * pPtr->wt();

            if (i == 1)
            {
                TEqvTargetCurve_[superCelli] = i;

                TEqvTmp = totT;
            }
            else if (totT < TEqvTmp)
            {
                TEqvTargetCurve_[superCelli] = i;

                TEqvTmp = totT;
            }
        }

        if (!independentYEqvCurve_)
        {
            YEqvTargetCurve_[superCelli] = TEqvTargetCurve_[superCelli];
        }
        else
        {
            //- Find equivalent species curve independent of equivalent temperature curve.
            interpolationLookUpTable<scalar> YEqvInterpTab
            (
                "YEqv" + primarySpecies_,
                this->owner().mesh().time().constant(),
                this->owner().mesh()
            );

            scalar YEqvTmp;

            label nsPrime = this->owner().composition().slgThermo().carrier().species()[primarySpecies_];

            for (label i = 1; i < YEqvInterpTab.output().size(); i++)
            {
                scalar totY = 0;

                for (particleType* pPtr : cellParticles)
                    totY += sqr
                        (
                            pPtr->Y()[nsPrime] 
                          - YEqvInterpTab.lookUp(pPtr->XiC(this->cVarName()))[i]
                        ) * pPtr->wt();

                if (i == 1)
                {
                    YEqvTargetCurve_[superCelli] = i;

                    YEqvTmp = totY;
                }
                else if(totY < YEqvTmp)
                {
                    YEqvTargetCurve_[superCelli] = i;

                    YEqvTmp = totY;
                }
            }
        }
    }
}


template<class CloudType>
void Foam::FlameletCurves<CloudType>::EqvETargetValues
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{    
    // Reference to the particle manager
    const particleNumberController& pManager = this->owner().pManager();
   
   
    //- Find target curves
    this->findEquivalentTempSpeciesCurves();

    //- Interpolate values on the target curves
    interpolationLookUpTable<scalar> TEqvInterpTab
    (
        "TEqv",
        this->owner().mesh().time().constant(),
        this->owner().mesh()
    );

    forAll(YEqvETarget,specieI)
    {
        if (!this->solveEqvSpecie()[specieI])
        {
            YEqvETarget[specieI].primitiveFieldRef() = 0;
        }
        else
        {
            interpolationLookUpTable<scalar> YEqvInterpTab
            (
                "YEqv" + YEqvETarget[specieI].name(),
                this->owner().mesh().time().constant(),
                this->owner().mesh()
            );

            forAll(this->owner().mesh().cells(),celli)
            {
                label superCelli = pManager.getSuperCellID(celli);
                
                //- Get the index of the coupling Variable in the full set 
                //- of mmc Variables
        
                const HashTable<label, word>& indexOfCVinXi = XiC.cVarInXi();
                const word CVariable = this->cVarName();
        
                scalar XiE = XiC.Vars(indexOfCVinXi[CVariable]).field()[celli];//!!!
                
                YEqvETarget[specieI][celli] = 
                YEqvInterpTab.lookUp(XiE)[YEqvTargetCurve_[superCelli]];
                
                TEqvETarget[celli] = 
                TEqvInterpTab.lookUp(XiE)[TEqvTargetCurve_[superCelli]];
            }
        }
    }
}

// ************************************************************************* //

