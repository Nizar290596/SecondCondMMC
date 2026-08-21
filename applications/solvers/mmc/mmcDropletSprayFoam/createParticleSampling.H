/*---------------------------------------------------------------------------*\
|                        Create Particle Statistics                           |
|-----------------------------------------------------------------------------|
| Creats a particle sampling object used for the particle statistics          |
| in post processing                                                          |
\*---------------------------------------------------------------------------*/

    IOdictionary cloudProperties
    (
        IOobject
        (
            "cloudProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    particleStatistics pStats
    (
        cloudProperties.subDict("particleSampling"),
        mesh
    );

    IOdictionary fuelCloudProperties
    (
        IOobject
        (
            "fuelCloudProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    particleStatistics pStatsFuelCloud
    (
        fuelCloudProperties.subDict("particleSampling"),
        mesh
    );


