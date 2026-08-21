// Create the fuel cloud

if (pSets.size() > 1)
    FatalError 
        << "mmcDropletSprayDNSFoam only works for a single mmc cloud!"
        << exit(FatalError);


basicDropletSprayDNSThermoCloud fuelCloud
(
    "fuelCloud",
    rhoFiltered,
    UFiltered,
    g,
    pSets[0].slgThermo(),
    true
);

fuelCloud.setMMCCloud(pSets[0]);

