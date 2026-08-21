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
 
Description
    Exact discrete-time integration of the Ornstein-Uhlenbeck (OU) process
    used by the second-conditioning extension of the MMC framework.
 
    The OU process (unit stationary variance):
 
        dω = -(1/τ)·ω dt  +  √(2/τ) dW
 
    Exact update (valid for any time step, no truncation error):
 
        ω(t+Δt) = ω(t)·exp(-Δt/τ)  +  √(1 − exp(-2Δt/τ)) · ξ
 
    where  ξ ~ N(0,1)  and  τ = tauOU.
 
    Stationary distribution: ω ~ N(0,1), independent of Δt.
 
    Two entry points are provided:
 
      - OUStateUpdate()        : primitive scalar update, no cloud/particle
                                 dependencies.
      - applyOUProcessUpdate() : template helper that checks secondCondFlag,
                                 draws the Gaussian variate from the cloud RNG,
                                 and writes the result back to the particle.
 
    Note: beta (cloud.secondCondBeta()) is NOT used here.  It enters only
    when computing the modified conditioning variable
        φ° = φ · exp(β · ω_OU)
    which is handled in a separate step once the base scalar φ is defined.
 
SourceFiles
    -  (header-only, all functions inline)
 
\*---------------------------------------------------------------------------*/
 
#ifndef OUProcessUpdate_H
#define OUProcessUpdate_H
 
#include "scalar.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// Primitive exact OU update — no cloud/particle dependencies.
//
// Parameters
//   omegaOld : current OU state ω(t)
//   dt       : time step Δt  [s]
//   tauOU    : OU relaxation time τ  [s]
//   xi       : N(0,1) variate for this step
//
// Returns ω(t+Δt).
inline scalar OUStateUpdate
(
    const scalar omegaOld,
    const scalar dt,
    const scalar tauOU,
    const scalar xi
)
{
    const scalar decay = exp(-dt / tauOU);
    return omegaOld * decay + sqrt(1.0 - decay*decay) * xi;
}
 
 
// Template helper: advance the OU state of a single particle by one time step.
//
// Called from MixingPopeParticle::setCellValues() (once per global Δt).
// Skips particles with secondCondFlag == 0 or tauOU <= SMALL.
template<class ParticleType, class CloudType>
inline void applyOUProcessUpdate
(
    ParticleType& p,
    CloudType& cloud,
    const scalar dt
)
{
    if (p.secondCondFlag() == 1)
    {
        const scalar tauOU = cloud.secondCondTauOU();
 
        if (tauOU > SMALL)
        {
            const scalar xi = cloud.rndGen().Normal(0, 1);
            p.omegaOU() = OUStateUpdate(p.omegaOU(), dt, tauOU, xi);
        }
    }
}
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //

