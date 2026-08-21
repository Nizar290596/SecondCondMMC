// Create the fuel cloud

// Only works for a single mmc cloud
if (pSets.size() > 1)
    FatalError 
        << "mmcDropletSprayFoam only works for a single mmc cloud!"
        << exit(FatalError);


basicDropletSprayThermoCloud fuelCloud
(
    "fuelCloud",
    rho,
    U,
    g,
    pSets[0].slgThermo(),
    true
);

// Set the mmcCloud
fuelCloud.setMMCCloud(pSets[0]);
