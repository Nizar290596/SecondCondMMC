    Info << nl << "Reading thermophysical properties" << endl;

    // Create the reaction thermo model
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    basicMultiComponentMixture& composition = thermo.composition();

    PtrList<volScalarField>& YEqvE = composition.Y();


    //- Equivalent fields for thermophysical coupling
    autoPtr<psiReactionThermo> pEqvThermo(psiReactionThermo::New(mesh));

    psiReactionThermo& thermoEqvETarget = pEqvThermo();

    PtrList<volScalarField>& YEqvETarget = thermoEqvETarget.composition().Y();

    volScalarField& TEqvETarget = thermoEqvETarget.T();

    volScalarField hEqvETarget = thermoEqvETarget.he();

    //- MMC variables (reference and coupling)

        //- Define location of file with mmcVariables definitions
        const fileName solverDir =
            mmcStatusMessage::solverBasePath("applications/solvers/mmc/mmcSprayFoam");

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

        mmcVarSet XiTarget(entries,mesh);


    // - Other quantities

    const word inertSpecie(thermo.lookup("inertSpecie"));

    Info<< "Creating field p\n" << endl;
    volScalarField& p = thermo.p();

    Info<< "Creating field hEqvE\n" << endl;
    volScalarField& hEqvE = thermo.he();

    Info<< "Creating field psi\n" << endl;
    const volScalarField& psi = thermo.psi();

    Info<< "Creating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //- Field to check density coupling targets
    volScalarField TTARGET
    (
        IOobject
        (
            "TTARGET",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        TEqvETarget
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

    Info<< "Creating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    Info<< "Creating field ghf\n" << endl;
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    #include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );




    IOdictionary transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar Sc_(readScalar(transportProperties_.lookup("SchmidtNumber")));
    Info << "schmidt number is set to be: " << Sc_ << endl;

    scalar Sct_(readScalar(transportProperties_.lookup("TurbulentSchmidtNumber")));
    Info << "turbulent schmidt number is set to be: " << Sct_ << endl;

    //- a_ISO
    Info<< "Creating field D for a_ISO\n" << endl;
    volScalarField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        (turbulence->nu()/Sc_)
    );

    //- a_ISO
    Info<< "Creating field Dt for a_ISO\n" << endl;
    volScalarField Dt
    (
        IOobject
        (
            "Dt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
           // IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        (turbulence->nut()/Sct_)
    );

    Info<< "Creating field DEff\n" << endl;
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
        D + Dt
    );


    Info << "Creating field mu" << nl << endl;
    volScalarField mu("mu", turbulence->mu());

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


    Info << "Creating field nu" << nl << endl;
    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->mu()/rho
    );

    Info << "Creating field nut" << nl << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->mut()/rho
    );

    Info << "Creating field epsilon" << nl << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (turbulence->epsilon())/rho*rho
    );

    Info << "Creating field k" << nl << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (turbulence->k())/rho*rho
    );


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

    volScalarField evaporationRate
    (
        IOobject
        (
            "evaporationRate",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 0, -1, 0, 0), 0)
    ); 

    // Calculating D10(AMD) and D32(SMD) [m]
    // *************************************
    Info<< "Creating Mean Diameters variables fields\n" << endl;

    // Arithmetic Mean Diamter (D10)
    volScalarField numAMD  //Sum(Ni*Di) = D10 Numerator
    (
        IOobject
        (
            "D10Num",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength, 0)
    );

    volScalarField denAMD //Sum(Ni) = D10 Denominator
    (
        IOobject
        (
            "D10Den",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    // Sauter Mean Diamter (D32)
    volScalarField numSMD //Sum(Ni*Di^3) = D32 Numerator
    (
        IOobject
        (
            "D32Num",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0), 0)
    );

    volScalarField denSMD //Sum(Ni*Di^2) = D32 Denominator
    (
        IOobject
        (
            "D32Den",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, 0, 0, 0), 0)
    );
    // -------------------------------------------------------------


    //- Mass Flux of the lagrangian particle [m^3/(m^2*s)]
    Info<< "Creating field particleFlux\n" << endl;
    volScalarField particleFlux
    (
        IOobject
        (
            "particleFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    );

    volScalarField particleFlux10
    (
        IOobject
        (
            "particleFlux10",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    );

    volScalarField particleFlux20
    (
        IOobject
        (
            "particleFlux20",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    );

    volScalarField particleFlux30
    (
        IOobject
        (
            "particleFlux30",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    );
        
    volScalarField particleFlux40
    (
        IOobject
        (
            "particleFlux40",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
    );

    volScalarField particleFlux50
    (
        IOobject
        (
            "particleFlux50",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity, 0)
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

    if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
    {
        const volScalarField& delta =
            mesh.objectRegistry::lookupObject<volScalarField>("delta");

       // DeltaE.internalField() = delta.internalField();
        DeltaE = delta;
    }


    Info<< "Creating field muEff\n" << endl;
    volScalarField muEff
    (
        IOobject
        (
            "muEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->muEff()
    );


    Info<< "Creating field gradRho\n" << endl;
    volVectorField gradRho
    (
        IOobject
        (
            "gradRho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(rho * DEff)
    );


    Info<< "Creating field magSqrGradf\n" << endl;
    volScalarField magSqrGradf
    (
        IOobject
        (
            "magSqrGradf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSqr(fvc::grad(Xi.Vars("f").field()))
    );

//    Info<< "Creating field boundaryRho\n" <<endl;
//    volScalarField& boundaryRho = rho;

    #include "createMRF.H"

    //- Addition to treat convection term in PDEs of YEqvE and hEqvE as
    //- well as for Xi
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsY;

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsXi;

    forAll(Xi, i)
    {
        fieldsXi.add(Xi.Vars(i).field());
    }

