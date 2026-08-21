/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "EnvelopeModel.H"
#include "mathematicalConstants.H"
#include "meshTools.H"
#include "volFields.H"

using namespace Foam;//::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EnvelopeModel<CloudType>::EnvelopeModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
//**//    f1_(this->owner().linkFG().size()+1),
//**//    f0_(this->owner().linkFG().size())
    f1_(), //**//
    f0_() //**//

{}
    


template<class CloudType>
Foam::EnvelopeModel<CloudType>::EnvelopeModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
//**//    f1_(this->owner().linkFG().size()+1),
//**//    f0_(this->owner().linkFG().size())
    f1_(), //**//
    f0_() //**//
{}

template<class CloudType>
Foam::EnvelopeModel<CloudType>::EnvelopeModel
(
    const EnvelopeModel<CloudType>& efm
)
:
    CloudSubModelBase<CloudType>(efm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EnvelopeModel<CloudType>::~EnvelopeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
//**//inline Foam::scalarField& Foam::EnvelopeModel<CloudType>::f1()
inline Foam::scalar& Foam::EnvelopeModel<CloudType>::f1() //**//
{
    return f1_;
}

template<class CloudType>
//**//inline Foam::scalarField& Foam::EnvelopeModel<CloudType>::f0()
inline Foam::scalar& Foam::EnvelopeModel<CloudType>::f0() //**//
{
    return f0_;
}

template<class CloudType> // global switches (per default turned off)
void Foam::EnvelopeModel<CloudType>::setSwitches()
{
    TIgnition_ = this->coeffDict().lookupOrDefault("TIgnition", 0.0);
    Info << "Ignition Temperature for envelope flames is set to " << TIgnition_ << "K" << endl;
}

template<class CloudType>
void Foam::EnvelopeModel<CloudType>::initialize
(
    const scalarField Y0,
    const scalar T0,
    const label sCell
)
{}

template<class CloudType>
bool Foam::EnvelopeModel<CloudType>::checkState
(
    const scalar T0
)
{
    //TODO: implement checks for slipe velocity and droplet size (see switches)
    
    bool envelopeSwitch = true;
    if(T0 < TIgnition_) envelopeSwitch = false;
    return envelopeSwitch;
}

template<class CloudType> // set boundary values for mixture fraction
void Foam::EnvelopeModel<CloudType>::setBoundaries
(
    const scalar f0,
    const scalar f1
)
{
    f0_ = f0;
    f1_ = f1;
}

template<class CloudType> // current implemntation sets f equal to sum(YXLiq)=1 (only fuel component in liquid)
//**//Foam::scalarField Foam::EnvelopeModel<CloudType>::fNoEnvelope
Foam::scalar Foam::EnvelopeModel<CloudType>::fNoEnvelope //**// can be extended that droplets carry their own fX information
(
    const scalarField Y0,
    const scalarField Y1,
    const scalarField YX
)
{
//**//    scalarField fX;
    scalarField linkFG = this->owner().linkFG();
    scalar linkFGSize = linkFG.size();
//**//    fX.setSize(linkFGSize+1,-1.0);
    scalar fX = -1.0; //**//
    scalar fSum = 0.0;
//**//    f1_[linkFGSize] = 0.0;
//**//    f1_ = 0.0; //**//

//**//    forAll(linkFG,nfs)
//**//    {
//**//        f1_[nfs] = Y1[nfs];          //individual liquid mixture fractions equal the corresponding Mass fractions
//**//        f1_[linkFGSize] += f1_[nfs]; //overall liquid mixture fraction is the sum over all individual mixture fractions
//**//        f1_ += Y1[nfs]; //**//
//**//    }

    forAll(linkFG,nfs)
    {
//**//        f0_[nfs] = f0 * f1_[nfs];              //calculate individual f0 by scaling the overall f0 by the individual f1
//**//        fX[nfs] = (f1_[nfs] - f0_[nfs]) / (Y1[nfs] - Y0[linkFG[nfs]] + VSMALL) * (YX[nfs] - Y0[linkFG[nfs]]) + f0_[nfs];
        fX = (f1_ - f0_) / (Y1[nfs] - Y0[linkFG[nfs]] + VSMALL) * (YX[nfs] - Y0[linkFG[nfs]]) + f0_; //**//
//**//        fSum += fX[nfs];
        fSum += fX; //**//
    }
//**//    fX[linkFGSize] = fSum;
    fX = fSum / linkFGSize; //**//

    return fX;
}

template<class CloudType>
Foam::scalarField Foam::EnvelopeModel<CloudType>::YNoEnvelope
(
    const scalarField Y0,
    const scalarField Y1,
    const scalar fX
) const
{
    scalarField YX;
    YX.setSize(this->owner().pope().composition().componentNames().size(),0.0);
//**//    scalar linkFGSize = this->owner().linkFG().size();
    scalar sumY = 0.0;
    forAll(this->owner().pope().composition().componentNames(), ns)
    {
//**//        YX[ns] = (Y1[ns] - Y0[ns]) / (f1_[linkFGSize] - f0 + VSMALL) * (fX - f0) + Y0[ns];
        YX[ns] = (Y1[ns] - Y0[ns]) / (f1_ - f0_ + VSMALL) * (fX - f0_) + Y0[ns]; //**//
        sumY += YX[ns];
    }
    YX /= sumY;

    return YX;
}

template<class CloudType>
Foam::scalar Foam::EnvelopeModel<CloudType>::haNoEnvelope
(
    const scalarField YX,
    const scalar p,
    const scalar TD
) const
{
    scalar haX;
    haX = this->owner().pope().composition().particleMixture(YX).Ha(p, TD);

    return haX;
}

template<class CloudType>
Foam::scalar Foam::EnvelopeModel<CloudType>::TNoEnvelope
(
     const scalar TD
) const
{
    return TD;
}

template<class CloudType>
//**//Foam::scalarField Foam::EnvelopeModel<CloudType>::fSurface
Foam::scalar Foam::EnvelopeModel<CloudType>::fSurface //**//
(
    const scalarField Y0,
    const scalarField Y1,
    const scalar f0,
    const scalarField YSurf
)
{
    notImplemented
    (
        "fSurface return function not implemented"
    );
//**//    scalarField fDummy;
//**//    fDummy.setSize(this->owner().linkFG().size(),0.0);
//**//    return fDummy;
    return -2.0; //**//
}

template<class CloudType>
Foam::scalarField Foam::EnvelopeModel<CloudType>::YSurface
(
    const scalarField Y0Liq,
    const scalarField Y1,
    const scalar fSurf
)
{
    notImplemented
    (
        "YSurface return function not implemented"
    );
    scalarField YDummy;
    YDummy.setSize(this->owner().linkFG().size(),0.0);

    return YDummy;
}

template<class CloudType>
Foam::scalarField Foam::EnvelopeModel<CloudType>::YFilm
(
    const scalarField Y0Liq,
    const scalarField Y1,
    const scalar fFilm
) const
{
    notImplemented
    (
        "YFilm return function not implemented"
    );
    scalarField YDummy;
    YDummy.setSize(this->owner().linkFG().size(),0.0);

    return YDummy;
}

template<class CloudType>
Foam::scalar Foam::EnvelopeModel<CloudType>::haSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar ha0,
    const scalar ha1,
    const scalar fSurf
) const
{
    notImplemented
    (
        "haSurface return function not implemented"
    );

    return 0;
}

template<class CloudType>
Foam::scalar Foam::EnvelopeModel<CloudType>::TSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar TGas,
    const scalar fSurf,
    const scalar haSurf
) const
{
    notImplemented
    (
        "TSurface return function not implemented"
    );

    return 0;
}

template<class CloudType> //- turns composition of n_fuel components into n_thermo components 
Foam::scalarField Foam::EnvelopeModel<CloudType>::Y1
(
    const scalarField YFuelVap
) const
{
    scalarField Y1;
    Y1.setSize(this->owner().pope().composition().componentNames().size(), 0.0);
    scalar sumY = 0.0;
    forAll(this->owner().linkFG(), nfs)
    {
        Y1[this->owner().linkFG()[nfs]] = YFuelVap[nfs];
        sumY += YFuelVap[nfs];
    }
    Y1 /= sumY;
    
    return Y1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EnvelopeModelNew.C"


// ************************************************************************* //

