/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "basicSprayThermoCloud.H"

#include "makeParcelCloudFunctionObjects.H"
 
// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeThermoParcelSurfaceFilmModels.H" // thermo variant
 
// Thermodynamic
#include "makeParcelHeatTransferModels.H"
 
// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// Special srpay models 
#include "makeSprayParcelInjectionModels.H"
#include "makeSprayParcelForces.H"
#include "makeParcelEvaporationModels.H"
#include "makeParcelCpModels.H"
#include "makeParcelEnvelopeModels.H"
#include "makeParcelFPRadiationModels.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Default OpenFOAM sub models

    makeParcelCloudFunctionObjects(basicSprayThermoCloud);

    // Kinematic sub-models
    makeThermoParcelForces(basicSprayThermoCloud);
    makeParcelDispersionModels(basicSprayThermoCloud);
    makeParcelInjectionModels(basicSprayThermoCloud);
    makeParcelPatchInteractionModels(basicSprayThermoCloud);
    makeParcelStochasticCollisionModels(basicSprayThermoCloud);
    makeParcelSurfaceFilmModels(basicSprayThermoCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicSprayThermoCloud);


    // MPPIC sub-models
    makeMPPICParcelDampingModels(basicSprayThermoCloud);
    makeMPPICParcelIsotropyModels(basicSprayThermoCloud);
    makeMPPICParcelPackingModels(basicSprayThermoCloud);


// Added submodels of spray model

    makeSprayParcelInjectionModels(basicSprayThermoCloud);
    makeSprayParcelForces(basicSprayThermoCloud);
    makeParcelCpModels(basicSprayThermoCloud);
    makeParcelEnvelopeModels(basicSprayThermoCloud);
    makeParcelEvaporationModels(basicSprayThermoCloud);
    makeParcelFPRadiationModels(basicSprayThermoCloud);


// ************************************************************************* //
