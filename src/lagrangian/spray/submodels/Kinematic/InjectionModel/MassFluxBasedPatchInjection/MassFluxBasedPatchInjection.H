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
    Foam::MassFluxBasedPatchInjection

Description
    Patch injection.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Initial parcel velocity
      - Injection volume flow rate

    Properties:
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    MassFluxBasedPatchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef MassFluxBasedPatchInjection_H
#define MassFluxBasedPatchInjection_H

#include "InjectionModel.H"
#include "massFluxBasedPatchInjectionBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class MassFluxBasedPatchInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MassFluxBasedPatchInjection
:
    public InjectionModel<CloudType>,
    public massFluxBasedPatchInjectionBase
{
    // Private data

        //- Injection duration [s]
        scalar duration_;
        label patchId_;
        //- Number of parcels to introduce per second []
//        const label parcelsPerSecond_;

        scalar averageField_;

        scalar dTRun_;

        scalar maxMassFlux_;

        label facei_;

//        fvPatchField<scalar> massFluxPatch_;


//        const Switch sprayAngle_;
//        const scalar maxAngle_;
        const scalar maxRadius_;

        const point spacialInjectionOffset_;

         // mass flux stuff

         const word fieldName_;

         tmp<volScalarField> massFlux_;
//         surfaceScalarField massFluxPatch_;

         const scalar averageParcelMass_;

         const scalar massFlowRateScalingFactor_;

         const Switch stokesScaling_;

         const scalar L0_;
         //- Initial parcel velocity [m/s]
        const scalar UGasMax_;

        const scalar ULiqMin_;

        const scalar StkMin_;

        const scalar StkMax_;

        //- Flow rate profile relative to SOI []
        const TimeFunction1<scalar> flowRateProfile_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("massFluxBasedPatchInjection");


    // Constructors

        //- Construct from dictionary
        MassFluxBasedPatchInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        MassFluxBasedPatchInjection(const MassFluxBasedPatchInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new MassFluxBasedPatchInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MassFluxBasedPatchInjection();


    // Member Functions

        //- Inherit updateMesh from massFluxBasedPatchInjectionBase
        using massFluxBasedPatchInjectionBase::updateMesh;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from massFluxBasedPatchInjectionBase
            using massFluxBasedPatchInjectionBase::setPositionAndCell;

            //- Set the injection position and owner cell, tetFace and tetPt
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

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel//,
//                const vector U
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
    #include "MassFluxBasedPatchInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
