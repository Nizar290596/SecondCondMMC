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

#include "ThermoPhysicalCouplingModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoPhysicalCouplingModel<CloudType>::setThermoPhysicalCoupling()
{
    solveEqvSpecie_.setSize(this->owner().composition().componentNames().size(), false);

    // Read in which species are to be used in the coupling, and the primary
    // species for dynamic matching
    const dictionary couplingDict(this->owner().cloudProperties().subDict("thermophysicalCoupling"));

    forAll(this->owner().composition().componentNames(), ns)
    {
        if (couplingDict.found(this->owner().composition().componentNames()[ns]))
        {
            solveEqvSpecie_[ns] = true;
        }
    }

    couplingDict.lookup("condVariable") >> cVarName_;

    // This entry is consumed as XiC(cVarName_) by ReactingPopeParticle::calc()
    // and by FlameletCurves - it is the coupling variable handed to the
    // reaction model as the mixture fraction. It must therefore name a
    // registered coupling variable. Check it here: otherwise the first offence
    // is a HashTable::at() abort part-way through particle tracking, on one
    // arbitrary rank, long after start-up.
    //
    // To condition a coupling model on something else (e.g. the modified
    // progress variable), set condVariable in that model's own coefficients
    // sub-dictionary instead - see KernelEstimation.
    if (!XiC_.cVarInXiC().found(cVarName_))
    {
        FatalErrorInFunction
            << "thermophysicalCoupling/condVariable is " << cVarName_
            << ", which is not a registered coupling variable." << nl
            << "Registered coupling variables: " << XiCNames_ << nl
            << "This entry names the coupling variable passed to the reaction "
            << "model as the mixture fraction, so it must be one of those."
            << exit(FatalError);
    }


    tauRelaxBlending_ = couplingDict.lookupOrDefault<Switch>("tauRelaxBlending",false);

    const scalar tauTmpTarget(readScalar(couplingDict.lookup("tauRelax")));

    tauRelaxTarget_.value() = tauTmpTarget;

    if (tauRelaxBlending_)
    {
        const scalar tauTmpStart(readScalar(couplingDict.lookup("tauRelaxStart")));

        tauRelaxStart_.value() = tauTmpStart;

        const scalar tauTmpDelta(readScalar(couplingDict.lookup("tauRelaxDelta")));

        tauRelaxDelta_ = tauTmpDelta;
    }
    else
    {
        tauRelaxStart_.value() = tauTmpTarget;
        tauRelaxDelta_ = 1;
    }

    const word tUnits(couplingDict.lookup("tauUnits"));

    tauUnits_ = tUnits;

    if (tauUnits_ != "time" && tauUnits_ != "timestep")
    {
        FatalErrorIn
        (
            "void Foam::ThermoPhysicalCouplingModel<CloudType>::setThermoPhysicalCoupling()"
        )
        << "Illegal value for keyword tauUnits" << nl
        << "Legal values are - time - or - timestep -" <<nl
        << abort(FatalError);
    }
    
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoPhysicalCouplingModel<CloudType>::ThermoPhysicalCouplingModel
(
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    SubModelBase<CloudType>(owner),

    XiC_(Xi),

    XiCNames_(Xi.cVarInXi().toc()),

    Indicator_(0.0*Xi.Vars(0).field())
{
    setThermoPhysicalCoupling();
}


template<class CloudType>
Foam::ThermoPhysicalCouplingModel<CloudType>::ThermoPhysicalCouplingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const mmcVarSet& Xi
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type),

    XiC_(Xi),

    XiCNames_(Xi.cVarInXi().toc()),

    Indicator_(0.0*Xi.Vars(0).field())
{
    setThermoPhysicalCoupling();
}


template<class CloudType>
Foam::ThermoPhysicalCouplingModel<CloudType>::ThermoPhysicalCouplingModel
(
    const ThermoPhysicalCouplingModel<CloudType>& cm
)
:
    SubModelBase<CloudType>(cm),

    XiC_(cm.XiC()),

    XiCNames_(cm.XiCNames()),

    Indicator_(cm.Indicator()),
    
    cVarName_(cm.cVarName()),
    
    solveEqvSpecie_(cm.solveEqvSpecie_),
        
    tauRelaxBlending_(cm.tauRelaxBlending_),
        
    tauRelaxStart_(cm.tauRelaxStart_),
        
    tauRelaxTarget_(cm.tauRelaxTarget_),
        
    tauRelaxDelta_(cm.tauRelaxDelta_),
        
    tauUnits_(cm.tauUnits_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoPhysicalCouplingModel<CloudType>::~ThermoPhysicalCouplingModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoPhysicalCouplingModelNew.C"

// ************************************************************************* //

