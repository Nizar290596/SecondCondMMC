/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "basicPremixedReactingPopeCloud.H"

// Kinematic
#include "makeItoPopeParticleInflowBoundaryModels.H"

// Thermodynamic
#include "makeThermoPopeParticleCompositionModels.H"
#include "makeThermoPopeParticleThermoPhysicalCouplingModels.H"

// Soot Model
//#include "makeSootPopeParticleSootModels.H"
//#include "makeSootPopeParticleRadiationModels.H"

// Aerosol                                                                                                               
//#include "makeAerosolPopeParticleSynthesisModels.H" 

// Mixing
#include "makePremixedMixingPopeParticleMixingModels.H"

// Reacting
#include "makeReactingPopeParticleReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    makeItoPopeParticleInflowBoundaryModels(basicPremixedReactingPopeCloud);

    // Thermo sub-models
    makeThermoPopeParticleCompositionModels(basicPremixedReactingPopeCloud);
    makeThermoPopeParticleThermoPhysicalCouplingModels(basicPremixedReactingPopeCloud);

    // Mixing sub-models
    makePremixedMixingPopeParticleMixingModels(basicPremixedReactingPopeCloud);

    // Reacting sub-models
    makeReactingPopeParticleReactionModels(basicPremixedReactingPopeCloud);
    
    // Soot sub-models
//    makeSootPopeParticleSootModels(basicPremixedReactingPopeCloud);
//    makeSootPopeParticleCloudRadiationModels(basicPremixedReactingPopeCloud);

}


// ************************************************************************* //
