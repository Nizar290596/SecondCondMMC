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

#include "MassFluxBasedPatchInjection.H"
#include "TimeFunction1.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MassFluxBasedPatchInjection<CloudType>::MassFluxBasedPatchInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    massFluxBasedPatchInjectionBase
    (
        owner.mesh(), 
        this->coeffDict().template get<word>("patchName"), 
        readScalar(this->coeffDict().lookup("maxRadius"))
    ),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    maxRadius_(this->coeffDict().lookupOrDefault("maxRadius",1.0)),
    spacialInjectionOffset_(this->coeffDict().lookup("spacialInjectionOffset")),
    fieldName_(this->coeffDict().lookup("fieldName")),
    massFlux_
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,
                owner.db().time().timeName(),
                owner.db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            owner.mesh()
        )
    ),
    averageParcelMass_(readScalar(this->coeffDict().lookup("averageParcelMass"))),
    massFlowRateScalingFactor_
    (
        this->coeffDict().lookupOrDefault("massFlowRateScalingFactor",1.0)
    ),
    stokesScaling_(this->coeffDict().lookup("stokesScaling")),
    L0_(this->coeffDict().lookupOrDefault("L0",1.0)),
    UGasMax_(this->coeffDict().lookupOrDefault("UGasMax",1.0)),
    ULiqMin_(this->coeffDict().lookupOrDefault("ULiqMin",1.0)),
    StkMin_(this->coeffDict().lookupOrDefault("StkMin",1.0)),
    StkMax_(this->coeffDict().lookupOrDefault("StkMax",1.0)),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    )
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    massFluxBasedPatchInjectionBase::updateMesh(owner.mesh());
    patchId_ = this->owner().mesh().boundaryMesh().findPatchID(patchName_);

    // Set total volume/mass to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);
    Info << "Stokes Scaling: " << stokesScaling_ << endl;
    if(stokesScaling_)
    {
        Info << "    Stk(min/max):    " << StkMin_ << "   " << StkMax_ << endl;
        Info << "    UGasMax/ULiqMin: " << UGasMax_ << "   " << ULiqMin_ << endl;
        Info << "    L0:              " << L0_ << endl;
    }
    Info << "massFlowRateScalingFactor: " << massFlowRateScalingFactor_ << endl;
}


template<class CloudType>
Foam::MassFluxBasedPatchInjection<CloudType>::MassFluxBasedPatchInjection
(
    const MassFluxBasedPatchInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    massFluxBasedPatchInjectionBase(im),
    duration_(im.duration_),
    maxRadius_(im.maxRadius_),
    spacialInjectionOffset_(im.spacialInjectionOffset_),
    fieldName_(im.fieldName_),
    massFlux_(im.massFlux_),
    averageParcelMass_(im.averageParcelMass_),
    massFlowRateScalingFactor_(im.massFlowRateScalingFactor_),
    stokesScaling_(im.stokesScaling_),
    L0_(im.L0_),
    UGasMax_(im.UGasMax_),
    ULiqMin_(im.ULiqMin_),
    StkMin_(im.StkMin_),
    StkMax_(im.StkMax_),
    flowRateProfile_(im.flowRateProfile_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MassFluxBasedPatchInjection<CloudType>::~MassFluxBasedPatchInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MassFluxBasedPatchInjection<CloudType>::updateMesh()
{
    massFluxBasedPatchInjectionBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::MassFluxBasedPatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::MassFluxBasedPatchInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    volScalarField& massFlux = massFlux_.ref();
    massFlux.correctBoundaryConditions();
    const polyPatch& cPatch = this->owner().mesh().boundaryMesh()[patchId_];

    const surfaceScalarField& magSf = this->owner().mesh().magSf();

//    massFluxPatch_=massFlux.boundaryField()[patchId];

    scalar patchArea = 0.0;
    scalar massFlowRate = 0.0;
    scalar deltaT = this->owner().mesh().time().deltaTValue();//time1 - time0;

    maxMassFlux_ = 0.0;

    forAll(cPatch, faceI)
    {
        scalar massFluxI = massFlux.boundaryField()[patchId_][faceI];
        scalar faceAreaI = magSf.boundaryField()[patchId_][faceI];
        patchArea += faceAreaI;
        massFlowRate += faceAreaI*massFluxI;
        if(massFluxI>maxMassFlux_)
            maxMassFlux_=massFluxI;
    }
    reduce(patchArea, sumOp<scalar>());
    reduce(massFlowRate, sumOp<scalar>());

    massFlowRate*=massFlowRateScalingFactor_;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar mass = deltaT*massFlowRate;

        scalar nParcels = mass/averageParcelMass_;

        Info << "deltaT: " << deltaT << ", mass: " << mass << endl;

        Random& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcels > 0
         && (
               nParcels - scalar(nParcelsToInject)
             > rnd.globalPosition(scalar(0), scalar(1))
            )
        )
        {
            ++nParcelsToInject;
        }

//        scalar massGlobal = mass;
//        reduce(massGlobal, sumOp<scalar>());
        averageField_ = averageField_*dTRun_ + mass;
        dTRun_ += deltaT;
        averageField_ /= dTRun_;

        Info << "Patch Area: " << patchArea << ", mass flow rate scaled / unscaled: " << massFlowRate << " / " << massFlowRate/massFlowRateScalingFactor_ << ", average: " << averageField_ << ", dTRun: " << dTRun_ << endl;
        Info << "mass flow rates (unscaled / scaled / average): " << massFlowRate/massFlowRateScalingFactor_ << "   " << massFlowRate << "   " << averageField_ << endl;
        Info << "massFlowRate and mass: " << this->owner().mesh().time().timeName() << "   " << massFlowRate << "   "  << mass << endl; 
	return nParcelsToInject;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::MassFluxBasedPatchInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::MassFluxBasedPatchInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    massFluxBasedPatchInjectionBase::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        massFlux_.ref(),
        maxMassFlux_,
        position,
        cellOwner,
        tetFacei,
        tetPti,
	facei_
    );

    position+=spacialInjectionOffset_;
}


template<class CloudType>
void Foam::MassFluxBasedPatchInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel)
{
    // set particle velocity
//    scalar vFraction = readScalar(this->coeffDict().lookup("localVelocityFraction"));


//    volVectorField& U = this->owner().U().ref();
//    const polyPatch& cPatch = this->owner().mesh().boundaryMesh()[patchId_];
//    parcel.U() = this->owner().U()[parcel.cell()] * vFraction + U0_ * (1-vFraction);
    // set particle diameter
    parcel.d() = sizeDistribution_->sample();

    if(!stokesScaling_)
        parcel.U() = this->owner().U().boundaryField()[patchId_][facei_];
    else
    {
        scalar muGas = this->owner().mu()[parcel.cell()];
        scalar tauP = parcel.rho()*pow(parcel.d(),2)/(18*muGas+VSMALL);
        scalar tauG = L0_/(mag(this->owner().U().boundaryField()[patchId_][facei_])+UGasMax_)*2.0;
        scalar stokes = tauP/tauG;
        vector UMax = this->owner().U().boundaryField()[patchId_][facei_];
        scalar reductionFactor = mag(this->owner().U().boundaryField()[patchId_][facei_])/ULiqMin_;
        vector UMin = this->owner().U().boundaryField()[patchId_][facei_] / reductionFactor;

        if(stokes<=StkMin_)
            parcel.U() = UMax;
        else if(stokes>=StkMax_)
            parcel.U() = UMin;
        else
        {
            parcel.U() = (UMin - UMax) / (StkMax_ - StkMin_) * (stokes - StkMin_) + UMax;
        }

        Pout << "D: " << parcel.d() 
             << ", UGas: " << this->owner().U().boundaryField()[patchId_][facei_]
             << ", Stk: " << stokes
             << ", Uini: " << parcel.U() << endl;
    }

/*
    scalar muGas = this->owner().mu()[parcel.cell()];
    Pout << "muGas: " << muGas << endl;
    scalar tauP = parcel.rho()*pow(parcel.d(),2)/(18*muGas+VSMALL);
    scalar tauG0 = characteristicL_/characteristicU_;
    scalar tauG1 = characteristicL_/(mag(this->owner().U().boundaryField()[patchId_][facei_])+VSMALL);
    scalar tauG2 = characteristicL_/(mag(this->owner().U().boundaryField()[patchId_][facei_])+characteristicU_)*2.0;
    scalar stokes0 = tauP/tauG0;
    scalar stokes1 = tauP/tauG1;
    scalar stokes2 = tauP/tauG2;
*/
//    parcel.U() = this->owner().U().boundaryField()[patchId_][facei_];//* vFraction_ + U0_ * (1-vFraction_);

/*
     scalar momentum = pow(parcel.d(),3)/6*3.14159*parcel.rho()*mag(parcel.U());
     Pout << "grepStokes: " << parcel.d() << token::TAB << stokes0 << token::TAB << stokes1 << token::TAB << stokes2 <<endl;   
     Pout << "momentum [kgm/s]: " << momentum << endl;
*/
    //    parcel.initCellFacePt();
}


template<class CloudType>
bool Foam::MassFluxBasedPatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::MassFluxBasedPatchInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
