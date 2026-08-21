// Create volVectorField U
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensioned<vector>("u",dimVelocity,vector(0,0,0))
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho0",dimDensity,1)
);

volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho0",dimless,1E-9)
);

volScalarField Dt
(
    IOobject
    (
        "Dt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho0",dimless,1E-9)
);


volScalarField DEff
(
    IOobject
    (
        "DEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho0",dimless,1E-9)
);


volScalarField DeltaE
(
    IOobject
    (
        "DeltaE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DeltaE", dimLength, SMALL)
);


volVectorField gradRhoDEff
(
    IOobject
    (
        "gradRhoDEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(rho * DEff)
);


//- Define location of file with mmcVariables definitions
const fileName solverDir = 
            mmcStatusMessage::solverBasePath("applications/solvers/mmc/mmcFoam"); 

const IOdictionary control
(
    IOobject
    (
        "mmcVariablesDefinitions",
        // Location of solver
        solverDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


const dictionary entries(control.subDict("variables"));



mmcVarSet Xi(entries,mesh);
