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
    Foam::KernelEstimation

Description
    Templated parcel single phase, multi-component class

SourceFiles
    KernelEstimation.C

\*---------------------------------------------------------------------------*/

#ifndef KernelEstimation_H
#define KernelEstimation_H

#include "ThermoPhysicalCouplingModel.H"
#include "kdTree.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef List<scalar> densParticle;

/*---------------------------------------------------------------------------*\
                    Class KernelEstimation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KernelEstimation
:
    public ThermoPhysicalCouplingModel<CloudType>
{
    //- comparison functor used for sorting in different directions
    class lessArg
    {
        label ii_;

        public:

            lessArg(label ii) {ii_=ii;}

            inline bool operator()
            (
                densParticle * one,
                densParticle * two
            ) const
            {
                return (*one)[ii_]< (*two)[ii_];
            };

            inline bool operator()
            (
                densParticle & one,
                densParticle & two
            ) const
            {
                return (one)[ii_]< (two)[ii_];
            };

            inline bool operator()
            (
                const scalar bound,
                densParticle * two
            ) const
            {
                return bound < (*two)[ii_];
            };

            inline bool operator()
            (
                densParticle * one,
                const scalar bound
            ) const
            {
                return (*one)[ii_] < bound;
            };
    };

    // Private data

        const fvMesh& mesh_;

        labelList Yindexes_;
        label numYEqv_;
        label numXiC_;

        const label nI_;
        const label nx_;
        const label ny_;
        const label nz_;
        const label nwt_;
        const label nT_;
        const label nDist_;
        const label nYs_;

        label nYe_;
        label nXiCs_;
        label nXiCe_;

        label numP_;

        const scalar fm_;
        const scalar rMaxMax_;

	//- Number of nearest neighbours used to build each cell kernel
        const label nNearest_;

        DynamicList<densParticle>  particleList_;
        DynamicList<densParticle>  LESList_;

        List<scalar> maxVal_;
        List<scalar> minVal_;

        List<scalar> maxXiCVal_;
        List<scalar> minXiCVal_;

        const scalar fLow_;
        const scalar fHigh_;

        const word inertSpecie_;

        const scalar dfMax_;
        const scalar C2_;

        Switch debug_;

	//- Condition on the modified progress variable phi-degree (a
        //  particle-only quantity) instead of a registered Eulerian coupling
        //  field. Enabled when condVariable == "phiModified".
        bool phiModEnabled_;

        //- Array slot of the carrier coupling variable whose position is reused
        //  to hold the conditioning value (in the XiC subset and full Xi set).
        label condSlotXiC_;
        label condSlotXi_;

        //- Per-cell phi-degree projected from the flagged particles; only
        //  built/used in phiModEnabled_ mode.
        autoPtr<volScalarField> phiModCell_;

        void buildParticleList();
        void buildLESParticleList();

	//- Project the flagged particles' phi-degree onto the mesh
        void buildPhiModCell();

        void computeTargets
        (
            const mmcVarSet& XiC,
            PtrList<volScalarField>& YEqvETarget,
            volScalarField& TEqvETarget
        );

//        void mergeLists();
public:

    //- Runtime type information
    TypeName("KernelEstimation");

    //- Type of particle the cloud was instantiated for
    typedef typename CloudType::particleType particleType;

    // Constructors

    //- Construct from dictionary
    KernelEstimation(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

    //- Destructor
    virtual ~KernelEstimation();

    // Member Functions

    // Access

    // Mixture properties

    //- Sources

        //- Equivalent temperature and species target values
        virtual void EqvETargetValues
        (
            const mmcVarSet& Xi,
            PtrList<volScalarField>& YEqvETarget,
            volScalarField& TEqvETarget
        );

  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "KernelEstimation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
