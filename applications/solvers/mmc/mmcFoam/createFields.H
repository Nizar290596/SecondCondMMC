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

        mmcVarSet XiTarget(entries,mesh);

    // - Other quantities

    const word inertSpecie(thermo.lookup("inertSpecie"));

    Info << "Creating field p" << nl << endl;
    volScalarField& p = thermo.p();

    Info << "Creating field hEqvE" << nl << endl;
    volScalarField& hEqvE = thermo.he();

    Info << "Creating field psi" << nl << endl;
    const volScalarField& psi = thermo.psi();

    Info << "Creating field rho" << nl << endl;
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

    Info << "Reading field U" << nl << endl;
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

    Info << "Creating field dpdt" << nl << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info << "Creating field kinetic energy K" << nl << endl;
    volScalarField K("K", 0.5*magSqr(U));

    Info << "Creating field g.h" << nl << endl;
    volScalarField gh("gh", g & mesh.C());

    Info << "Creating field ghf" << nl << endl;
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    #include "compressibleCreatePhi.H"

    Info << "Creating turbulence model" << nl << endl;
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
    Info << "Schmidt number is set to be: " << Sc_ << endl;

    scalar Sct_(readScalar(transportProperties_.lookup("TurbulentSchmidtNumber")));
    Info << "Turbulent schmidt number is set to be: " << Sct_ << nl << endl;

    //- aISO
    Info << "Creating field D for aISO" << nl << endl;
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

    //- aISO
    Info << "Creating field Dt for aISO" << nl << endl;
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

    Info << "Creating field DEff" << nl << endl;
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


    Info << "Creating field gradRho" << nl << endl;
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

    #include "createMRF.H"

    //- Addition to treat convection term in PDEs of YEqvE and hEqvE as
    //- well as for Xi
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsY;

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsXi;

    forAll(Xi, i)
    {
        fieldsXi.add(Xi.Vars(i).field());
    }

