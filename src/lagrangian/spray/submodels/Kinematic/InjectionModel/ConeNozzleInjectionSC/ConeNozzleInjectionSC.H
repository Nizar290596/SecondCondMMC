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

Class
    Foam::ConeNozzleInjectionSC

Description
    Cone injection disc model with Size conditioned velocities input.

    User specifies:
      - time of start of injection
      - injector position
      - direction (along injection axis)
      - parcel flow rate
      - inner and outer half-cone angles

    Properties:
      - Parcel diameters obtained by size distribution model.

      - Parcel velocity is calculated as:
        - Constant velocity:
          \verbatim
          U = \<specified by user\>
          \endverbatim

SourceFiles
    ConeNozzleInjectionSC.C

\*---------------------------------------------------------------------------*/

#ifndef ConeNozzleInjectionSC_H
#define ConeNozzleInjectionSC_H

#include "InjectionModel.H"
#include "Vector.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class TimeFunction1;

class distributionModel;

/*---------------------------------------------------------------------------*\
                     Class ConeNozzleInjectionSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConeNozzleInjectionSC
:
    public InjectionModel<CloudType>
{
public:

    //- Injection method enumeration
    enum injectionMethod
    {
        imPoint,
        imDisc
    };

    //- Flow type enumeration
    enum flowType
    {
        ftConstantVelocity
    };


private:

    // Private data
         typedef VectorSpace<Vector<scalar>, scalar, 7> septuple;

        //- Point/disc injection method
        injectionMethod injectionMethod_;

        //- Flow type
        flowType flowType_;

        //- Outer nozzle diameter [m]
        const scalar outerDiameter_;

        //- Inner nozzle diameter [m]
        const scalar innerDiameter_;

        //- Injection duration [s]
        scalar duration_;

        //- Injector position [m]
        vector position_;

        //- Cell containing injector position []
        label injectorCell_;

        //- Index of tet face for injector cell
        label tetFacei_;

        //- Index of tet point for injector cell
        label tetPti_;

        //- Injector direction []
        vector direction_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeFunction1<scalar> flowRateProfile_;

        //- Inner half-cone angle relative to SOI [deg]
        const TimeFunction1<scalar> thetaInner_;

        //- Outer half-cone angle relative to SOI [deg]
        const TimeFunction1<scalar> thetaOuter_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> sizeDistribution_;

		// Add Random Fluctuations to the Input Velocities 
        // U = UMean + rand*URMS  = UMean (if fluctuations_ = FALSE)
        const Switch randFluctuations_;
		
		// Use Transverse Velocities of the Carrier phase rather
		// than using the experimental table inputs	
        const Switch cpTransVel_;


        // Tangential vectors to the direction vector

            //- First tangential vector
            vector tanVec1_;

            //- Second tangential vector
            vector tanVec2_;

            //- Injection vector orthogonal to direction
            vector normal_;


        // Velocity model coefficients

            //- Input Table
            List<septuple> UList_;

    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();

        //- Set the injection flow type
        void setFlowType();


public:

    //- Runtime type information
    TypeName("coneNozzleInjectionSC");


    // Constructors

        //- Construct from dictionary
        ConeNozzleInjectionSC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ConeNozzleInjectionSC(const ConeNozzleInjectionSC<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ConeNozzleInjectionSC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConeNozzleInjectionSC();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);

        //- Linear Interpolate
        inline scalar interpLinear(const scalar x,const scalar y1,const scalar y0,const scalar x1,const scalar x0)
        {
            const scalar m = (y1 - y0)/(x1 - x0);
            const scalar c = y0 - m*x0;
            return m*x + c;
        }
        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ConeNozzleInjectionSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
