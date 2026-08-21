    Info << nl << "Reading thermophysical properties" << endl;

    // Create the reaction thermo model
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    basicMultiComponentMixture& composition = thermo.composition();

    PtrList<volScalarField>& YEqvE = composition.Y();

    //Info << nl << "Thermophysical backward coupling is NOT implemented! "<< nl << endl;
    //- Equivalent fields for thermophysical coupling
    autoPtr<psiReactionThermo> pEqvThermo(psiReactionThermo::New(mesh));
    
    psiReactionThermo& thermoEqvETarget = pEqvThermo();

    PtrList<volScalarField>& YEqvETarget = thermoEqvETarget.composition().Y();

    volScalarField& TEqvETarget = thermoEqvETarget.T();

    volScalarField hEqvETarget = thermoEqvETarget.he(); 

    

    // - Other quantities

    const word inertSpecie(thermo.lookup("inertSpecie"));

    Info << "The inert Specie is" << inertSpecie << nl << endl;

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
    
    Info << "Creating field gammafit" << nl << endl;
    volScalarField gammafit
    (
        IOobject
        (
            "gammafit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gammafit", dimless, 1.0)
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
    Info << "Creating molecular diffusivity field D for aISO" << nl << endl;
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
    Info << "Creating turbulent diffusivity field Dt for aISO" << nl << endl;
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

    Info << "Creating field DEff (D+Dt)" << nl << endl;
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

	Info<< "Creating field lambda\n" << endl;
    volScalarField vb
    (
        IOobject
        (
            "vb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("vb", dimless, 1.0)
    );


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

//- MMC variables (reference and coupling)

        //- Define location of file with mmcVariables definitions
        //const fileName solverDir = 
        //    mmcStatusMessage::solverBasePath("applications/solvers/mmc/vlShadowPositionFoam"); 
        
        //const fileName userDir = getEnv("WM_PROJECT_USER_DIR");
        //const fileName mmcVersion = "mmcFoam";
        //const fileName solverDir = userDir/mmcVersion/"applications/solvers/mmc/vlShadowPositionFoam";

//	const fileName solverDir = "/home/562/nh4682/OpenFOAM/nh4682-v2406/mmcFoam/applications/solvers/mmc/vlShadowPositionFoam";
//	const fileName solverDir = "/home/nizar/OpenFOAM/nizar-v2406/mmcFoam/applications/solvers/mmc/SPFoam";
	const fileName solverDir = "/home/562/nh4682/OpenFOAM/nh4682-v2406/mmcFoam/applications/solvers/mmc/SPFoam";

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
        
        Info << "The mmcVariable definitions are:\n" << entries.toc() << endl;

        mmcVarSet Xi(entries,mesh);

        mmcVarSet XiTarget(entries,mesh);



    #include "createMRF.H"

    //- Addition to treat convection term in PDEs of YEqvE and hEqvE as
    //- well as for Xi
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsY;

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsXi;

    forAll(Xi, i)
    {
        fieldsXi.add(Xi.Vars(i).field());
    }
    
    
    Info << "Creating field uPrime" << nl << endl;
    volScalarField uPrime
    (
        IOobject
        (
            "uPrime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("uPrime", dimVelocity, 0.0)
    );
    
    //- Create uPrimek based on sgs turbulent kinetic energy
    Info<< "Creating field uPrimek\n" << endl;
    volScalarField uPrimek
    (
        IOobject
        (
            "uPrimek",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("uPrimek",dimVelocity, 0.0)
    );
    
    volScalarField dynsl0
    (
        IOobject
        (
            "dynsl0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dynsl0", dimVelocity, 2.5E-01)
    );
    
    volScalarField dyndeltal0
    (
        IOobject
        (
            "dyndeltal0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dyndeltal0", dimLength, 5.0E-04)
    );
    
    
    //- Create LESdelta to compute thickening factor
    Info << "Creating field LES_delta" << nl << endl;
    volScalarField LES_delta
    (
        IOobject
        (
            "LES_delta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("LES_delta", dimLength, 0.0005)
    );

    if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
    {
        const volScalarField& delta =
            mesh.objectRegistry::lookupObject<volScalarField>("delta");

        LES_delta.ref() = delta.internalField();
    }

    volVectorField laplaceSim
    (
        IOobject
        (
            "laplaceSim",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("laplaceSim", dimensionSet(0,-1,-1,0,0,0,0), vector::zero)
     );

    volVectorField Uint
    (
        IOobject
        (
            "Uint",
            runTime.timeName(),
            evalLaplaceMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        evalLaplaceMesh,
        dimensionedVector("Uint", dimVelocity, vector::zero)
    );

    volVectorField coarseLap
    (
        IOobject
        (
            "coarseLap",
            runTime.timeName(),
            evalLaplaceMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        evalLaplaceMesh,
        dimensionedVector("coarseLap",dimensionSet(0,-1,-1,0,0,0,0), vector::zero)
    );
    
        //- Create local Karlovitz number
    Info<< "Creating local Karlovitz number\n" << endl;
    volScalarField Ka
    (
        IOobject
        (
            "Ka",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("Ka", dimless, 1.0)
    );

    //- Create local Damkoehler number
    Info<< "Creating local Damkoehler number\n" << endl;
    volScalarField Da
    (
        IOobject
        (
            "Da",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("Da", dimless, 1.0)
    );
    
    //- Create local RL
    Info<< "Creating local RL \n" << endl;
    volScalarField RL
    (
        IOobject
        (
            "RL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("RL", dimless, 1.0)
    );

    //- Create local RV
    Info<< "Creating local RV \n" << endl;
    volScalarField RV
    (
        IOobject
        (
            "RV",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("RV", dimless, 1.0)
    );

    //- Create local a_vb
    Info<< "Creating local a_vb \n" << endl;
    volScalarField a_vb
    (
        IOobject
        (
            "a_vb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("a_vb", dimless, 1.0)
    );

    //- Create local fu
    Info<< "Creating local fu \n" << endl;
    volScalarField fu
    (
        IOobject
        (
            "fu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	//1.0
        dimensionedScalar("fu", dimless, 1.0)
    );

    //- Create local fdelta
    Info<< "Creating local fdelta \n" << endl;
    volScalarField fdelta
    (
        IOobject
        (
            "fdelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("fdelta", dimless, 1.0)
    );
    
    //- Create local Resgs
    Info<< "Creating local Resgs \n" << endl;
    volScalarField Resgs
    (
        IOobject
        (
            "Resgs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Resgs", dimless, 1.0)
    );   
    
    //- Create local fre
    Info<< "Creating local fre \n" << endl;
    volScalarField fre
    (
        IOobject
        (
            "fre",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("fre", dimless, 1.0)
    );
