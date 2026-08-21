/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
inline bool
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
variableTimeStep() const
{
    return variableTimeStep_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::label
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::timeSteps() const
{
    return timeSteps_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::logFile(const word& name) const
{
    mkDir(this->mesh().time().path()/"TDAC");
    return autoPtr<OFstream>
    (
        new OFstream
        (
            this->mesh().time().path()/"TDAC"/name
        )
    );
}


template<class ReactionThermo, class ThermoType>
inline Foam::PtrList<Foam::volScalarField>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::Y()
{
    return this->Y_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::autoPtr<Foam::particleChemistryReductionMethod<ReactionThermo, ThermoType>>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::mechRed()
{
    return mechRed_;
}


template<class ReactionThermo, class ThermoType>
inline void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::setActive
(
    const label i
)
{
    this->thermo().composition().setActive(i);
}


template<class ReactionThermo, class ThermoType>
inline bool Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::active
(
    const label i
) const
{
    return this->thermo().composition().active(i);
}


template<class ReactionThermo, class ThermoType>
inline void
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
setNsDAC(const label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class ReactionThermo, class ThermoType>
inline void
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}


template<class ReactionThermo, class ThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::Field<Foam::label>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class ThermoType>
inline const Foam::Field<Foam::label>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::Field<bool>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::reactionsDisabled()
{
    return reactionsDisabled_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalarField&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::completeC()
{
    return completeC_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalarField&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::List<Foam::List<Foam::specieElement>>&
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::specieComp()
{
    return specieComp_;
}


template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::resetTabulationResults()
{
    tabulationResults_ = 2;
}


template<class ReactionThermo, class ThermoType>
inline void
Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::
incrsTimeStepCounter()
{
    timeSteps_++;
}

template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::setZeroCPUTimeData()
{
    // CPU time analysis
    reduceMechCpuTime_ = 0;
    addNewLeafCpuTime_ = 0;
    growCpuTime_ = 0;
    solveChemistryCpuTime_ = 0;
    searchISATCpuTime_ = 0;

    // Average number of active species
    nActiveSpecies = 0;
    nAvg = 0;
}

template<class ReactionThermo, class ThermoType>
void Foam::particleTDACChemistryModel<ReactionThermo, ThermoType>::writeCPUTimeLog()
{

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        //Write the ISAT test log
        if(tabulation_->TestSwitch())
        tabulation_->ISATLifeStatistics();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (mechRed_->active() && nAvg && mechRed_->log())
    {
        // Write average number of species
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies/nAvg << endl;
    }

}





// ************************************************************************* //
