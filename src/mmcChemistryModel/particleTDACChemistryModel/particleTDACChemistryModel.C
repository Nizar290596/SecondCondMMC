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

#include "particleTDACChemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
particleTDACChemistryModel
(
	ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
    variableTimeStep_
    (
        this->mesh().time().controlDict().lookupOrDefault("adjustTimeStep", false)
     || fv::localEulerDdt::enabled(this->mesh())
    ),
    timeSteps_(0),
    dTlimitor_(this->lookupOrDefault("dTlimitor",1.0e8)),
    NsDAC_(this->nSpecie_),
    completeC_(this->nSpecie_, 0),
    reactionsDisabled_(this->reactions_.size(), false),
    specieComp_(this->nSpecie_),
    completeToSimplifiedIndex_(this->nSpecie_, -1),
    simplifiedToCompleteIndex_(this->nSpecie_),
    tabulationResults_(0)
{
    Info << "dTlimitor is " << dTlimitor_ << endl;
    
    basicMultiComponentMixture& composition = this->thermo().composition();

    // Store the species composition according to the species index
    speciesTable speciesTab = composition.species();

    // autoPtr to the specieCompositionTable
    auto specCompPtr = thermo.specieComposition();
    const auto specComp = specCompPtr();

    forAll(specieComp_, i)
    {
        specieComp_[i] = specComp[this->Y()[i].name()];
    }

    mechRed_ = particleChemistryReductionMethod<ReactionThermo, ThermoType>::New
    (
        *this,
        *this
    );

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialized (by default 'active' is true)
    if (mechRed_->active())
    {
        forAll(this->Y(), i)
        {
            IOobject header
            (
                this->Y()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i);
                this->Y()[i].writeOpt() = IOobject::NO_WRITE;
            }
        }
    }

    tabulation_ = particleChemistryTabulationMethod<ReactionThermo, ThermoType>::New
    (
        *this,
        *this
    );

    if (mechRed_->log())
    {
        cpuReduceFile_ = logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = logFile("nActiveSpecies.out");
    }

    if (tabulation_->log())
    {
        cpuAddFile_ = logFile("cpu_add.out");
        cpuGrowFile_ = logFile("cpu_grow.out");
        cpuRetrieveFile_ = logFile("cpu_retrieve.out");
    }


    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_ = logFile("cpu_solveCloud.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
~particleTDACChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(this->reactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType>& R = this->reactions_[i];

            scalar omegai = omega
            (
                R, c, T, p, pf, cf, lRef, pr, cr, rRef
            );

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl*omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr*omegai;
            }
        }
    }
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        // When using DAC, the ODE solver submit a reduced set of species the
        // complete set is used and only the species in the simplified mechanism
        // are updated
        this->c_ = completeC_;

        // Update the concentration of the species in the simplified mechanism
        // the other species remain the same and are used only for third-body
        // efficiencies
        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        for (label i=0; i<this->nSpecie(); i++)
        {
            this->c_[i] = max(0.0, c[i]);
        }
    }

    omega(this->c_, T, p, dcdt);

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W*this->c_[i];
    }

    scalar cp = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        // cp function returns [J/(kmol K)]
        cp += this->c_[i]*this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    // When mechanism reduction is active
    // dT is computed on the reduced set since dcdt is null
    // for species not involved in the simplified mechanism
    scalar dT = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        label si;
        if (reduced)
        {
            si = simplifiedToCompleteIndex_[i];
        }
        else
        {
            si = i;
        }

        // ha function returns [J/kmol]
        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    // limit the time-derivative, this is more stable for the ODE
    // solver when calculating the allowed time step
    const scalar dTLimited = min(dTlimitor_, mag(dT));
    dcdt[this->nSpecie_] = -dT*dTLimited/(mag(dT) + 1.0e-10);
    //dcdt[this->nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[this->nSpecie_ + 1] = 0;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarSquareMatrix& dfdc
) const
{
    const bool reduced = mechRed_->active();

    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the informations of the complete set
    // (i.e. for the third-body efficiencies)

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        this->c_ = completeC_;
        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        forAll(this->c_, i)
        {
            this->c_[i] = max(c[i], 0.0);
        }
    }

    dfdc = Zero;

    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions_[ri];

            const scalar kf0 = R.kf(p, T, this->c_);
            const scalar kr0 = R.kr(kf0, p, T, this->c_);

            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                if (reduced)
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kf = kf0;
                forAll(R.lhs(), i)
                {
                    const label si = R.lhs()[i].index;
                    const scalar el = R.lhs()[i].exponent;
                    if (i == j)
                    {
                        if (el < 1)
                        {
                            if (this->c_[si] > SMALL)
                            {
                                kf *= el*pow(this->c_[si] + VSMALL, el - 1);
                            }
                            else
                            {
                                kf = 0;
                            }
                        }
                        else
                        {
                            kf *= el*pow(this->c_[si], el - 1);
                        }
                    }
                    else
                    {
                        kf *= pow(this->c_[si], el);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (reduced)
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc(si, sj) -= sl*kf;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (reduced)
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc(si, sj) += sr*kf;
                }
            }

            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                if (reduced)
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kr = kr0;
                forAll(R.rhs(), i)
                {
                    const label si = R.rhs()[i].index;
                    const scalar er = R.rhs()[i].exponent;
                    if (i == j)
                    {
                        if (er < 1)
                        {
                            if (this->c_[si] > SMALL)
                            {
                                kr *= er*pow(this->c_[si] + VSMALL, er - 1);
                            }
                            else
                            {
                                kr = 0;
                            }
                        }
                        else
                        {
                            kr *= er*pow(this->c_[si], er - 1);
                        }
                    }
                    else
                    {
                        kr *= pow(this->c_[si], er);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (reduced)
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc(si, sj) += sl*kr;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (reduced)
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc(si, sj) -= sr*kr;
                }
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1e-3;

    omega(this->c_, T + delta, p, this->dcdt_);
    for (label i=0; i<this->nSpecie_; i++)
    {
        dfdc(i, this->nSpecie_) = this->dcdt_[i];
    }

    omega(this->c_, T - delta, p, this->dcdt_);
    for (label i=0; i<this->nSpecie_; i++)
    {
        dfdc(i, this->nSpecie_) =
            0.5*(dfdc(i, this->nSpecie_) - this->dcdt_[i])/delta;
    }

    dfdc(this->nSpecie_, this->nSpecie_) = 0;
    dfdc(this->nSpecie_ + 1, this->nSpecie_) = 0;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    jacobian(t, c, dfdc);

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    // Note: Uses the c_ field initialized by the call to jacobian above
    omega(this->c_, T, p, dcdt);
}

template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    notImplemented
        (
            "particleTDACChemistryModel::solve"
            "("
            "DeltaTType& "
            ")"
        );
    return (0);

}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    notImplemented
        (
            "particleTDACChemistryModel::solve"
            "("
            "scalar "
            ")"
        );
    return (0);
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    notImplemented
        (
            "particleTDACChemistryModel::solve"
            "("
            "scalarField& "
            ")"
        );
    return (0);
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    notImplemented
        (
            "particleTDACChemistryModel::solve"
            "("
            "scalarField&, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar"
            ") const"
        );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::particleCalculate
(
    const scalar t0,
    const scalar deltaT,
    const scalar hi,
    scalar pi,
    scalar Ti,
    scalarField& Y
)
{

    // Increment counter of time-step
    //timeSteps_++; //incremented in reactingCloud

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 1 : 0);

    //basicMultiComponentMixture& composition = this->thermo().composition();

    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();

    this->resetTabulationResults(); // FV field of tabulation actions -retrieve -grow -add. useless in this class

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    // Density of the mixture
    scalarField n(this->nSpecie_, 0.0);
    
    for (label i=0; i<this->nSpecie_; i++)
    {
        n[i] = Y[i]/this->specieThermo_[i].W();
    }
    
//    const scalar nTot = sum(n);
    
//    GasThermo mixture(0.0*this->specieThermo_[0]);
//    
//    for (label i=0; i<this->nSpecie_; i++)
//    {
//        mixture += (n[i]/nTot)*this->specieThermo_[i];
//    }

    ThermoType mixture
    (
        Y[0]*this->specieThermo_[0]
    );

    for (label i=1; i<this->nSpecie_; i++)
    {
        mixture += Y[i]*this->specieThermo_[i];
    }
    
    scalar rhoi=mixture.rho(pi,Ti);

//    scalarField c(this->nSpecie_, 0.0);
//    for (label i=0; i<this->nSpecie_; i++)
//    {
//        c[i] = rhoi*Y[i]/this->specieThermo_[i].W();
//    }

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    scalarField Rphiq(this->nEqns() + nAdditionalEqn);

    scalarField DIRphiq(this->nEqns() + nAdditionalEqn);

    for (label i=0; i<this->nSpecie_; i++)
    {
        c[i] = rhoi*Y[i]/this->specieThermo_[i].W();
        c0[i] = c[i]; //used to calculate RR
        phiq[i] = Y[i];
    }
    phiq[this->nSpecie()]=Ti;
    phiq[this->nSpecie() + 1]=pi;
    if (tabulation_->variableTimeStep())
    {
        phiq[this->nSpecie() + 2] = deltaT;
    }


    // Initialise time progress
    scalar timeLeft = deltaT;

    // Not sure if this is necessary
    Rphiq = Zero;

    clockTime_.timeIncrement();

    const label celli = 1;

    // When tabulation is active (short-circuit evaluation for retrieve)
    // It first tries to retrieve the solution of the system with the
    // information stored through the tabulation method
    if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
    {

        if(tabulation_->TestSwitch()) //calculate DI result to compare and get local error
        {
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    // completeC_ used in the overridden ODE methods
                    // to update only the active species
                    completeC_ = c;

                    // Solve the reduced set of ODE
                    this->solve
                    (
                        simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]
                    );

                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }
            for (label i=0; i<this->nSpecie_; i++)
            {
                DIRphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
            }
        }

        // Retrieved solution stored in Rphiq
        for (label i=0; i<this->nSpecie(); i++)
        {
            c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
        }

        searchISATCpuTime_ += clockTime_.timeIncrement();
    }
    // This position is reached when tabulation is not used OR
    // if the solution is not retrieved.
    // In the latter case, it adds the information to the tabulation
    // (it will either expand the current data or add a new stored point).
    else
    {
        // Store total time waiting to attribute to add or grow
        scalar timeTmp = clockTime_.timeIncrement();

        if (reduced)
        {
            // Reduce mechanism change the number of species (only active)
            mechRed_->reduceMechanism(c, Ti, pi);
            nActiveSpecies += mechRed_->NsSimp(); //should accumulate 
            nAvg++; //should be how many times accumulated
            scalar timeIncr = clockTime_.timeIncrement();
            reduceMechCpuTime_ += timeIncr;
            timeTmp += timeIncr;
        }

        // Calculate the chemical source terms
        while (timeLeft > SMALL)
        {
            scalar dt = timeLeft;
            if (reduced)
            {
                // completeC_ used in the overridden ODE methods
                // to update only the active species
                completeC_ = c;

                // Solve the reduced set of ODE
                this->solve
                (
                    simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]
                );

                for (label i=0; i<NsDAC_; i++)
                {
                    c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                }
            }
            else
            {
                this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
            }
            timeLeft -= dt;
        }

        {
            scalar timeIncr = clockTime_.timeIncrement();
            solveChemistryCpuTime_ += timeIncr;
            timeTmp += timeIncr;

        }

        // If tabulation is used, we add the information computed here to
        // the stored points (either expand or add)
        if (tabulation_->active())
        {
            forAll(c, i)
            {
                Rphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
            }
            if (tabulation_->variableTimeStep())
            {
                Rphiq[Rphiq.size()-3] = Ti;
                Rphiq[Rphiq.size()-2] = pi;
                Rphiq[Rphiq.size()-1] = deltaT;
            }
            else
            {
                Rphiq[Rphiq.size()-2] = Ti;
                Rphiq[Rphiq.size()-1] = pi;
            }
            label growOrAdd =
                tabulation_->add(phiq, Rphiq, rhoi, deltaT);

            if (growOrAdd)
            {
                this->setTabulationResultsAdd();
                addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
            else
            {
                this->setTabulationResultsGrow();
                growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
            }
        }

        // When operations are done and if mechanism reduction is active,
        // the number of species (which also affects nEqns) is set back
        // to the total number of species (stored in the mechRed object)
        if (reduced)
        {
            this->nSpecie_ = mechRed_->nSpecie();
        }
        deltaTMin = min(this->deltaTChem_[1], deltaTMin);
    }

    // Update particle mass fractions and temperature
//    const scalar cTot = sum(c);
//    
//    GasThermo mixture(0.0*this->specieThermo_[0]);

    for (label i=0; i<this->nSpecie_; i++)
    {
//        mixture += (c[i]/cTot)*this->specieThermo_[i];
        Y[i] = c[i]*this->specieThermo_[i].W()/rhoi;
    }


    return deltaTMin;
}

template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsAdd()
{
    tabulationResults_ = 0.0;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsGrow()
{
    tabulationResults_ = 1.0;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsRetrieve()
{
    tabulationResults_ = 2.0;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar 
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::getCPUTimeReaction()
{

    scalar CPUtime = solveChemistryCpuTime_;

    if(tabulation_->active())
    {
        CPUtime = searchISATCpuTime_ + growCpuTime_ + addNewLeafCpuTime_;
        
    }

    if(mechRed_->active())
    {
        CPUtime = CPUtime + reduceMechCpuTime_;
    }

    return CPUtime;

}



// ************************************************************************* //
