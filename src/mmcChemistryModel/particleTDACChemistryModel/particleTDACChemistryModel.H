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

Class
    Foam::particleTDACChemistryModel
 
Description
    Extends the TDACChemistryModel of OpenFOAM with a particleCalculate
    function that solves the chemistry on a particle level. 

SourceFiles
    particleTDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef particleTDACChemistryModel_H
#define particleTDACChemistryModel_H

#include "StandardChemistryModel.H"
#include "particleChemistryReductionMethod.H"
#include "particleChemistryTabulationMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class particleTDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class particleTDACChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    // Private member data

        bool variableTimeStep_;

        label timeSteps_;

        scalar dTlimitor_;

    // CPU time analysis
        scalar reduceMechCpuTime_ = 0;
        scalar addNewLeafCpuTime_ = 0;
        scalar growCpuTime_ = 0;
        scalar solveChemistryCpuTime_ = 0;
        scalar searchISATCpuTime_ = 0;

    // Average number of active species
        scalar nActiveSpecies = 0;
        scalar nAvg = 0;

        // Mechanism reduction
        label NsDAC_;
        scalarField completeC_;
        scalarField simplifiedC_;
        Field<bool> reactionsDisabled_;
        List<List<specieElement>> specieComp_;
        Field<label> completeToSimplifiedIndex_;
        DynamicList<label> simplifiedToCompleteIndex_;
        autoPtr<particleChemistryReductionMethod<ReactionThermo, ThermoType>> mechRed_;

        // Tabulation
        autoPtr<particleChemistryTabulationMethod<ReactionThermo, ThermoType>> tabulation_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;

        //- Log file for the average time spent adding tabulated data
        autoPtr<OFstream> cpuAddFile_;

        //- Log file for the average time spent growing tabulated data
        autoPtr<OFstream> cpuGrowFile_;

        //- Log file for the average time spent retrieving tabulated data
        autoPtr<OFstream> cpuRetrieveFile_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        //- Variable containing information about tabulation:
        //  0 -> add (direct integration)
        //  1 -> grow
        //  2 -> retrieve
        label tabulationResults_;


    // Private Member Functions

        //- Disallow copy constructor
        particleTDACChemistryModel(const particleTDACChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const particleTDACChemistryModel&);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("particleTDACChemistryModel");


    // Constructors

        //- Construct from mesh
        particleTDACChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~particleTDACChemistryModel();


    // Member Functions

        //- Return true if the time-step is variable and/or non-uniform
        inline bool variableTimeStep() const;

        //- Return the number of chemistry evaluations (used by ISAT)
        inline label timeSteps() const;

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        inline PtrList<volScalarField>& Y();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalarField& dcdt
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculate the estimation of chemical integration time step
        virtual scalar particleCalculate
        (
            const scalar t0,        // initial time
            const scalar deltaT,    // time step
            const scalar hi,        // enthalpy - currently not used
            scalar pi,              // pressure of particle i
            scalar Ti,              // temperature of particle i
            scalarField& Y          // list of species of particle i
        );

        // Chemistry model functions (overriding functions in
        // chemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);


        // ODE functions (overriding functions in chemistryModel to take into
        // account the variable number of species)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            //- Pure jacobian function for tabulation
            void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void solve
            (
                scalarField& c,
                scalar& T,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;


        // Mechanism reduction access functions

            inline void setNsDAC(const label newNsDAC);

            inline void setNSpecie(const label newNs);

            inline scalarField& completeC();

            inline scalarField& simplifiedC();

            inline Field<bool>& reactionsDisabled();

            inline bool active(const label i) const;

            inline void setActive(const label i);

            inline DynamicList<label>& simplifiedToCompleteIndex();

            inline Field<label>& completeToSimplifiedIndex();

            inline const Field<label>& completeToSimplifiedIndex() const;

            inline List<List<specieElement>>& specieComp();

            inline autoPtr<particleChemistryReductionMethod<ReactionThermo, ThermoType>>&
                mechRed();

            void setTabulationResultsAdd();

            void setTabulationResultsGrow();

            void setTabulationResultsRetrieve();

            inline void resetTabulationResults();

            inline void incrsTimeStepCounter();

            inline void setZeroCPUTimeData();

            inline void writeCPUTimeLog();

            scalar getCPUTimeReaction();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleTDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleTDACChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
