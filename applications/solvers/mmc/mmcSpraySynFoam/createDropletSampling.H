/*---------------------------------------------------------------------------*\
|                        Create Particle Statistics                           |
|-----------------------------------------------------------------------------|
| Creats a particle sampling object used for the fuel particle statistics     |
| in post processing                                                          |
| Note: This uses a different dictionary than the particle statistics         |
|       object.                                                               |
\*---------------------------------------------------------------------------*/

    IOdictionary fuelProperties
    (
        IOobject
        (
            "fuelCloud1Properties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    particleStatistics dStats
    (
        fuelProperties.subDict("dropletSampling"),
        mesh
    );
