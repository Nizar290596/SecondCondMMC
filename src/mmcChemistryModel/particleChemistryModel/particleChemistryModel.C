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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/
#include "particleChemistryModel.H"

template<class ReactionThermo, class ThermoType>
Foam::particleChemistryModel<ReactionThermo, ThermoType>::particleChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
{}


template<class ReactionThermo, class ThermoType>
Foam::scalar
Foam::particleChemistryModel<ReactionThermo, ThermoType>::particleCalculate
(
    const scalar t0,
    const scalar deltaT,
    const scalar hi,
    scalar pi,
    scalar Ti,
    scalarField& Y
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;
    
    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    if (Ti <= this->Treact_)
    {
        return deltaTMin;
    }

    // Store latest estimation of integration time step at
    // the cell location celli
    const label celli = 1;
    
    // Density of the mixture

    ThermoType mixture
    (
        Y[0]*this->specieThermo_[0]
    );

    for (label i=1; i<this->nSpecie_; i++)
    {
        mixture += Y[i]*this->specieThermo_[i];
    }

    scalar rhoi=mixture.rho(pi,Ti);

//    scalarField c(nSpecie_, 0.0);
//    scalarField c0(nSpecie_);

    for (label i=0; i<this->nSpecie_; i++)
    {
        this->c_[i] = rhoi*Y[i]/this->specieThermo_[i].W();
//        c0[i] = c_[i];
    }

    // Initialise time progress
    scalar timeLeft = deltaT;

    // Calculate the chemical source terms
    while (timeLeft > SMALL)
    {
        scalar dt = timeLeft;
        this->solve(this->c_, Ti, pi, dt, this->deltaTChem_[celli]);
        timeLeft -= dt;

        // Update the temperature

        rhoi= 0.0;

        for (label i=0; i<this->nSpecie_; i++)
        {
            rhoi += this->c_[i]*this->specieThermo_[i].W();
        }        

        for (label i=0; i<this->nSpecie_; i++)
        {
            Y[i] = this->c_[i]*this->specieThermo_[i].W()/rhoi;
        }
        
        ThermoType mixture
        (
            Y[0]*this->specieThermo_[0]
        );

        for (label i=1; i<this->nSpecie_; i++)
        {
            mixture += Y[i]*this->specieThermo_[i];
        }
/*
        ThermoType mixture
        (
            (specieThermo_[0].W()*c_[0])*specieThermo_[0]
        );

        for (label i=1; i<nSpecie_; i++)
        {
            mixture += (specieThermo_[i].W()*c_[i])*specieThermo_[i];
        }
        
        Ti = mixture.THa(hi, pi, Ti);
*/
    }


/*
    rhoi=mixture.rho(pi,Ti);

    for (label i=0; i<nSpecie_; i++)
    {
        Y[i] = c_[i]*specieThermo_[i].W()/rhoi;
    }
*/
    deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

    return deltaTMin;
}

