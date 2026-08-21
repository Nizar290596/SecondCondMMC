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
    Foam::TwoEquationModel

Description
    Dummy class for 'none' option - will raise an error if any functions are
    called that require return values.

SourceFiles
    TwoEquationModel.C

\*---------------------------------------------------------------------------*/

#ifndef TwoEquationSoot_H
#define TwoEquationSoot_H

#include "SootModel.H"
#include "ThermoPopeParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class TwoEquationModel Declaration
\*---------------------------------------------------------------------------*/
template<class CloudType>
class TwoEquationModel
:
    public SootModel<CloudType>
{
    //- Involved gas phase species
    //  by default, C2H2, O2, OH
    wordList species_;

    //- source for ySoot
    scalar dySootdt_;

    //- source for nSoot
    scalar dnSootdt_;

    //- soot diameter
    scalar diaSoot_;

    //- soot surface area
    scalar surfSoot_;

    //- soot kinetics coeffs includes:
    //  ----------------------------------------------
    //  0. inception        1. cogulation
    //  2. surface growth   3. oxidation
    //  Generally, for each list, the coeffs follows:
    //  [0] A, ie, the pre-factor
    //  [1] pressureDependency
    //  [2] TCoeff
    //  [3] TAct
    //  [4] SootSurfaceCoeff
    //  [5] CarbonNum
    //  ----------------------------------------------
    List<Tuple2<string, scalarList> > kineticsCoeffs_;

    //- soot density
    scalar rhoSoot_;

    //- scaling factor for nSoot_
    scalar nSootScaling_;

    //- number of carbon for inception
    label inceptionCarbon_;

    //- Switch for the Y source
    Switch gasSource_;

    //- gas source table
    //  1st column: reactionName_species
    //  2nd column: reaction coefficient
    List<Tuple2<string, scalar> > gasSourceTable_;

    List<Tuple2<string, scalar> > YsourceTable_;

    //- threshold for the soot intermittency
    scalar sootLimit_;

    //- Switch for turbulent mixing
    const Switch turbMixing_;

    //- debug mode
    Switch debug_;

    //- list of soot contributions
    //  ->[0] ySootInC2H2
    //  ->[1] ySootSgC2H2
    //  ->[2] ySootOxO2
    //  ->[3] ySootOxOH
    //  ->[4] nSootInC2H2
    //  ->[5] nSootCg
    scalarList sootContribs_;

public:

    //- Runtime type information
    TypeName("twoEquationModel");

    //- Constructors

        //- Construct from dictionary
        TwoEquationModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        TwoEquationModel(const TwoEquationModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SootModel<CloudType> > clone() const
        {
            return autoPtr<SootModel<CloudType> >
            (
                new TwoEquationModel<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~TwoEquationModel();

    //- Access the member data

        inline const wordList& species() const;

        inline const scalar& rhoSoot() const;

        inline const scalar& nSootScaling() const;

        inline const scalar& surfSoot() const;

        inline const List<Tuple2<string, scalarList> >& kineticsCoeffs() const;

        inline const label& inceptionCarbon() const;

        inline const Switch& gasSource() const;

        inline const List<Tuple2<string, scalar> >& gasSourceTable() const;

        inline const List<Tuple2<string, scalar> >& YsourceTable() const;

        inline const scalar& sootLimit() const;

        virtual Switch turbMixing() const override;

        inline const scalarList& sootContribs() const;

    //- Change the member data

        inline scalar& diaSoot();

        inline scalar& surfSoot();

        inline scalarList& sootContribs();

    //- species check, called before calculating soot volume fraction
    //  For two-equation soot model, normally (without considering gas sink)
    //  C2H2, O2 and OH are considered.
    //  When considering the gas phase sink/source,
    //  CO2, CO, H, and H2 then join in.
    void speciesCheck();

    //- Calculate the rate coeff. of cogulation
    scalar calRateCoeff
    (
        const scalar& kineticsCoeff,
        const scalar& T,
        const scalar& ySoot,
        const scalar& nSoot,
        const scalar& rho
    );

    //- Calculate the rate coeff. for inception, surface growth, and oxidation
    scalar calRateCoeff
    (
        const scalarList& kineticsCoeffs,
        const scalar& T
    );

    // Use to fill dySootdt_ and dnSootdt_
    // and by default
    // > inception species: C2H2
    // > surface growth species: C2H2
    // > oxidation, O2 and OH
    void fillInSource
    (
        const string& reactspecies,
        const scalarField& Y,
        const scalar& T,
        const scalar& ySoot,
        const scalar& nSoot,
        const scalar& rho,
        scalar& source4ySoot,
        scalar& source4nSoot
    );

    void updateSootVariables
    (
        const scalarField& Y,
        const scalar& T,
        scalar& ySoot,
        scalar& nSoot,
        const scalar& rho
    );

    template<class DataType>
    int indexInTable
    (
        const string& reactName,
        const List<Tuple2<string, DataType> >& Table
    );

    void fillInGasSourceTable
    (
        const string& reactName,
        const string& speciesTag,
        const scalar& W,
        const scalar& source,
        const scalar& deltaT,
        const scalar& rho
    );

    // update Ys
    void updateYs(scalarField& Y);

    //- updates soot variables: ySoot, nSoot, sootVf, sootIm
    //  and, transfer soot conContribs to particle
    //  ##### this function is called by each particle #####
    virtual void calSoot
    (
        scalarField& Y,
        const scalar& T,
        scalar& ySoot,
        scalar& nSoot,
        scalar& sootVf,
        scalar& sootIm,
        const scalar& rho,
        scalarList& sootContribs
    );

    //  listing the involved gas species for each step, does not work with soot calculations
    void speciesInfo();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TwoEquationModelI.H"

#ifdef NoRepository
#   include "TwoEquationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
