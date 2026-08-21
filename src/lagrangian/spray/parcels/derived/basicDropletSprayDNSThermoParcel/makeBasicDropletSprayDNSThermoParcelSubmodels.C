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

#include "basicDropletSprayDNSThermoCloud.H"

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

    makeParcelCloudFunctionObjects(basicDropletSprayDNSThermoCloud);

    // Kinematic sub-models
    makeThermoParcelForces(basicDropletSprayDNSThermoCloud);
    makeParcelDispersionModels(basicDropletSprayDNSThermoCloud);
    makeParcelInjectionModels(basicDropletSprayDNSThermoCloud);
    makeParcelPatchInteractionModels(basicDropletSprayDNSThermoCloud);
    makeParcelStochasticCollisionModels(basicDropletSprayDNSThermoCloud);
    makeParcelSurfaceFilmModels(basicDropletSprayDNSThermoCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicDropletSprayDNSThermoCloud);


    // MPPIC sub-models
    makeMPPICParcelDampingModels(basicDropletSprayDNSThermoCloud);
    makeMPPICParcelIsotropyModels(basicDropletSprayDNSThermoCloud);
    makeMPPICParcelPackingModels(basicDropletSprayDNSThermoCloud);

    // MMC coupling sub-models
    makeMMCCouplingModels(basicDropletSprayDNSThermoCloud, basicReactingPopeCloud);

    // Saturation pressure models
    makeLiquidPropertiesModels(basicDropletSprayDNSThermoCloud);

    // Added submodels of dropletspray model
    makeDropletSprayParcelInjectionModels(basicDropletSprayDNSThermoCloud);

// ************************************************************************* //
