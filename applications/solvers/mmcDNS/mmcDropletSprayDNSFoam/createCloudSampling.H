/*---------------------------------------------------------------------------*\
|                        Create Particle Statistics                           |
|-----------------------------------------------------------------------------|
| Creats a particle sampling object used for the particle statistics          |
| in post processing                                                          |
\*---------------------------------------------------------------------------*/

particleStatistics pStats
(
    pSets[0].cloudProperties(),
    mmcCloudMeshRef.cref()
);


particleStatistics pStatsFuelCloud
(
    fuelCloud.particleProperties().subDict("particleSampling"),
    mesh
);

