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

    // Resolve the registered coupling variable to use wherever an actual XiC
    // entry is required. condVariable may name a quantity that has no entry in
    // the mmcVarSet -- 'phiModified' is particle-only -- but the particle
    // chemistry and the flamelet look-ups still need a real coupling value, and
    // KernelEstimation needs a carrier array slot. Fall back to the first
    // registered coupling variable in that case.
    cVarNameXiC_ = word::null;

    forAll(XiCNames_, i)
    {
        if (XiCNames_[i] == cVarName_)
        {
            cVarNameXiC_ = cVarName_;
            break;
        }
    }

    if (cVarNameXiC_ == word::null)
    {
        if (XiCNames_.empty())
        {
            FatalErrorInFunction
                << "condVariable '" << cVarName_ << "' is not a registered "
                << "coupling variable, and no coupling variables are declared "
                << "in mmcVariablesDefinitions." << nl
                << "At least one couplingVar is required to carry the coupling "
                << "value for the particle chemistry."
                << exit(FatalError);
        }

        cVarNameXiC_ = XiCNames_[0];

        Info<< "ThermoPhysicalCouplingModel: condVariable '" << cVarName_
            << "' is not a registered coupling variable; using '"
            << cVarNameXiC_ << "' wherever an XiC entry is required "
            << "(particle chemistry, flamelet look-up, carrier slot)." << endl;
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

    cVarNameXiC_(cm.cVarNameXiC()),

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

