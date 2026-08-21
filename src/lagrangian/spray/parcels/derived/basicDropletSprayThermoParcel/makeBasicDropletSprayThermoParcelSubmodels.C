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

#include "basicDropletSprayThermoCloud.H"

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

// MMC coupling models
#include "makeMMCCouplingModels.H"
#include "basicReactingPopeCloud.H"

// Saturation pressure models
#include "makeLiquidPropertiesModels.H"

// Special dropletsrpay models 
#include "makeDropletSprayParcelInjectionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Default OpenFOAM sub models

    makeParcelCloudFunctionObjects(basicDropletSprayThermoCloud);

    // Kinematic sub-models
    makeThermoParcelForces(basicDropletSprayThermoCloud);
    makeParcelDispersionModels(basicDropletSprayThermoCloud);
    makeParcelInjectionModels(basicDropletSprayThermoCloud);
    makeParcelPatchInteractionModels(basicDropletSprayThermoCloud);
    makeParcelStochasticCollisionModels(basicDropletSprayThermoCloud);
    makeParcelSurfaceFilmModels(basicDropletSprayThermoCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicDropletSprayThermoCloud);


    // MPPIC sub-models
    makeMPPICParcelDampingModels(basicDropletSprayThermoCloud);
    makeMPPICParcelIsotropyModels(basicDropletSprayThermoCloud);
    makeMPPICParcelPackingModels(basicDropletSprayThermoCloud);

    // MMC coupling sub-models
    makeMMCCouplingModels(basicDropletSprayThermoCloud, basicReactingPopeCloud);

    // Saturation pressure models
    makeLiquidPropertiesModels(basicDropletSprayThermoCloud);

    // Added submodels of dropletspray model
    makeDropletSprayParcelInjectionModels(basicDropletSprayThermoCloud);

// ************************************************************************* //
