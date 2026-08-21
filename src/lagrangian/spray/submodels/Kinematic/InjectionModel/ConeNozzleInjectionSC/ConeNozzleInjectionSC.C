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

#include "ConeNozzleInjectionSC.H"
#include "TimeFunction1.H"
#include "mathematicalConstants.H"
#include "distributionModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjectionSC<CloudType>::setInjectionMethod()
{
    word injectionMethodType = this->coeffDict().getWord("injectionMethod");
    if (injectionMethodType == "disc")
    {
        injectionMethod_ = imDisc;
    }
    else if (injectionMethodType == "point")
    {
        injectionMethod_ = imPoint;

        // Set/cache the injector cell
        this->findCellAtPosition
        (
            injectorCell_,
            tetFacei_,
            tetPti_,
            position_,
            false
        );
    }
    else
    {
        FatalErrorInFunction
            << "injectionMethod must be either 'point' or 'disc'"
            << exit(FatalError);
    }
}

template<class CloudType>
void Foam::ConeNozzleInjectionSC<CloudType>::setFlowType()
{
    word flowType = this->coeffDict().getWord("flowType");
    if (flowType == "constantVelocity")
    {

        flowType_ = ftConstantVelocity;
    }

    else
    {
        FatalErrorInFunction
            << "flowType must be 'constantVelocity'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjectionSC<CloudType>::ConeNozzleInjectionSC
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_(imDisc),
    flowType_(ftConstantVelocity),
    outerDiameter_(readScalar(this->coeffDict().lookup("outerDiameter"))),
    innerDiameter_(readScalar(this->coeffDict().lookup("innerDiameter"))),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(this->coeffDict().lookup("position")),
    injectorCell_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    direction_(this->coeffDict().lookup("direction")),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "thetaOuter",
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
    ),
    randFluctuations_(this->coeffDict().lookupOrDefault("randomFluctuations",false)),
    cpTransVel_(this->coeffDict().lookupOrDefault("cpRadTanVel",false)),
    tanVec1_(Zero),
    tanVec2_(Zero),
    normal_(Zero),
    UList_(this->coeffDict().lookup("UList"))
{
    if (innerDiameter_ >= outerDiameter_)
    {
        FatalErrorInFunction
         << "innerNozzleDiameter >= outerNozzleDiameter" << nl
         << exit(FatalError);
    }
    Info << "Particle Injection Fluctuations: " << randFluctuations_ << endl;
    Info << "Take transverse velocity from gas phase: " << cpTransVel_ << endl;
    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    setFlowType();

    Random& rndGen = this->owner().rndGen();

    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    vector tangent = Zero;
    scalar magTangent = 0.0;

    while(magTangent < SMALL)
    {
        vector v = rndGen.sample01<vector>();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);

    updateMesh();
}


template<class CloudType>
Foam::ConeNozzleInjectionSC<CloudType>::ConeNozzleInjectionSC
(
    const ConeNozzleInjectionSC<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    flowType_(im.flowType_),
    outerDiameter_(im.outerDiameter_),
    innerDiameter_(im.innerDiameter_),
    duration_(im.duration_),
    position_(im.position_),
    injectorCell_(im.injectorCell_),
    tetFacei_(im.tetFacei_),
    tetPti_(im.tetPti_),
    direction_(im.direction_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    thetaInner_(im.thetaInner_),
    thetaOuter_(im.thetaOuter_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    randFluctuations_(im.randFluctuations_),
    cpTransVel_(im.cpTransVel_),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_),
    normal_(im.normal_),
    UList_(im.UList_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeNozzleInjectionSC<CloudType>::~ConeNozzleInjectionSC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConeNozzleInjectionSC<CloudType>::updateMesh()
{
    // Set/cache the injector cells
    switch (injectionMethod_)
    {
        case imPoint:
        {
            this->findCellAtPosition
            (
                injectorCell_,
                tetFacei_,
                tetPti_,
                position_
            );
        }
        default:
        {}
    }
}


template<class CloudType>
Foam::scalar Foam::ConeNozzleInjectionSC<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::ConeNozzleInjectionSC<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return floor((time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeNozzleInjectionSC<CloudType>::volumeToInject
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
void Foam::ConeNozzleInjectionSC<CloudType>::setPositionAndCell
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
    Random& rndGen = this->owner().rndGen();

    scalar beta = mathematical::twoPi*rndGen.sample01<scalar>();
    normal_ = tanVec1_*cos(beta) + tanVec2_*sin(beta);

    switch (injectionMethod_)
    {
        case imPoint:
        {
            position = position_;
            cellOwner = injectorCell_;
            tetFacei = tetFacei_;
            tetPti = tetPti_;

            break;
        }
        case imDisc:
        {
            scalar frac = rndGen.sample01<scalar>();
            scalar dr = outerDiameter_ - innerDiameter_;
            scalar r = 0.5*(innerDiameter_ + frac*dr);
            position = position_ + r*normal_;

            this->findCellAtPosition
            (
                cellOwner,
                tetFacei,
                tetPti,
                position,
                false
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
             << "Unknown injectionMethod type" << nl
             << exit(FatalError);
        }
    }
}


template<class CloudType>
void Foam::ConeNozzleInjectionSC<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    // OpenFOAM Random Library 
    Random& rndGen = this->owner().rndGen();

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();

    // Parcel Injection Postion (Radial)
    scalar pRadPos  = sqrt(pow(parcel.position().x(),2) + pow(parcel.position().y(),2));
    
    switch (flowType_)
    {
        case ftConstantVelocity:
        {
            // Find the value from the input table
            label n=0;
            for (label i=0; i<UList_.size(); i++)
            {
                if(parcel.d() > UList_[i][1] && parcel.d() <= UList_[i][2])
                {
                    label k = i;
                    while (UList_[k][0] <= pRadPos)
                    {
                        k++;
                    }
                    n=k; // Chosen value index ahead(++) of the interpolation loc)
                    break;
                }
            }

            if (n == 0)
            {
                FatalErrorInFunction
                     << "Parcel Sampled Diameter = " << parcel.d() << " is not in the input list." << nl
                     << exit(FatalError);
            }

            if (UList_[n][0] < pRadPos) //[TODO]: Not sure if this should exist or not.
            {
            FatalErrorInFunction
                 << "Parcel Position(n+1) = " << pRadPos << " does not have upper limit for radial location." << nl
                 << "r(Chosen) = " << UList_[n][0] << " and Table Index(Chosen) = " << nl
                 << "Check Outer Diamter of the injection model and radial limit of the table." << nl << nl
                 << exit(FatalError);
            }

            //Info << "n                      = " << n <<nl;
            //Info << "UMean(n)               = " << UList_[n][3] << nl;
            //Info << "UMean(n-1)             = " << UList_[n-1][3] << nl;
            //Info << "Parcel Diameter        = " << parcel.d() << nl;
            //Info << "Diameter Min   (Chosen)= " << UList_[n][1] << nl;
            //Info << "Diameter Max   (Chosen)= " << UList_[n][2] << nl;
            //Info << "RadLoc (n-1)           = " << UList_[n-1][0] << nl;
            //Info << "ParcelPos              = " << pRadPos << nl;    
            //Info << "RadLoc (n)             = " << UList_[n][0] << nl;
            
            // Random Fluctuations Components (U = UMean + U')
            scalar flucX = 0.0;
            scalar flucY = 0.0;
            scalar flucZ = 0.0;

            // Radial locations to interpolate velocities b/w
            const scalar x1 = UList_[n][0];
            const scalar x0 = UList_[n-1][0];

            // Interpolate the Mean and RMS values from the input table
            scalar UmeanAx  = interpLinear(pRadPos, UList_[n][3], UList_[n-1][3], x1, x0);
            scalar UmeanRad = interpLinear(pRadPos, UList_[n][4], UList_[n-1][4], x1, x0);
            scalar UrmsAx   = interpLinear(pRadPos, UList_[n][5], UList_[n-1][5], x1, x0);
            scalar UrmsRad  = interpLinear(pRadPos, UList_[n][6], UList_[n-1][6], x1, x0);

            //Info << "UmeanAx          = " << UmeanAx << nl << nl;

            if(randFluctuations_)
            {
                flucX = rndGen.GaussNormal<scalar>() * UrmsRad / sqrt(2.0);// * radius / parcel.position().x();
                flucY = rndGen.GaussNormal<scalar>() * UrmsRad / sqrt(2.0);// * radius / parcel.position().y();
                flucZ = rndGen.GaussNormal<scalar>() * UrmsAx;
            }

            if(cpTransVel_)
            {
				// Use carrier phase velocity for transverse velocities
                parcel.U() = this->owner().U()[parcel.cell()];
            }
            else
            {
				// Use table input values for transverse velocities
                parcel.U().x() = UmeanRad * parcel.position().x()/pRadPos + flucX;
                parcel.U().y() = UmeanRad * parcel.position().y()/pRadPos + flucY;
            }

            parcel.U().z() = UmeanAx + flucZ;            
            break;
        }

        default:
        {
        }
    }
}


template<class CloudType>
bool Foam::ConeNozzleInjectionSC<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ConeNozzleInjectionSC<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
