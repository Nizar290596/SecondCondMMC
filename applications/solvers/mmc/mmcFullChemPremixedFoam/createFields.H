    Info << nl << "Reading thermophysical properties" << endl;

    // Create the reaction thermo model
/*    autoPtr<BasicChemistryModel<psiReactionThermo>> pChemistry
    (
        BasicChemistryModel<psiReactionThermo>::New(thermo)
    );
    
    BasicChemistryModel<psiReactionThermo>& chemistry = pChemistry();

    psiReactionThermo& thermo = chemistry.thermo();*/
    
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    basicMultiComponentMixture& composition = thermo.composition();

    PtrList<volScalarField>& YEqvE = composition.Y();

    Info << nl << "Thermophysical backward coupling is NOT implemented! "<< nl << endl;
    //- Equivalent fields for thermophysical coupling
/*    autoPtr<psiReactionThermo> pEqvThermo(psiReactionThermo::New(mesh));
    
    psiReactionThermo& thermoEqvETarget = pEqvThermo();

    PtrList<volScalarField>& YEqvETarget = thermoEqvETarget.composition().Y();

    volScalarField& TEqvETarget = thermoEqvETarget.T();

    volScalarField hEqvETarget = thermoEqvETarget.he(); */

    //- MMC variables (reference and coupling)

        //- Define location of file with mmcVariables definitions
        const fileName solverDir = 
            mmcStatusMessage::solverBasePath("applications/solvers/mmc/mmcFullChemPremixedFoam"); 

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

//        mmcVarSet XiTarget(entries,mesh);

    // - Other quantities

    const word inertSpecie(thermo.lookup("inertSpecie"));

    Info << "Creating field p" << nl << endl;
    volScalarField& p = thermo.p();

    Info << "Creating field hEqvE" << nl << endl;
    volScalarField& hEqvE = thermo.he();
    
    Info << "Creating field T" << nl << endl;
    volScalarField& T = thermo.T();

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
/*    volScalarField TTARGET
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
    );*/

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
    
    forAll(YEqvE, i)
    {
        fieldsY.add(YEqvE[i]);
    }
    fieldsY.add(hEqvE);

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsXi;

    forAll(Xi, i)
    {
        fieldsXi.add(Xi.Vars(i).field());
    }
    
 // Fields for ATF
    
    //- Create flame sensor in every cell
    Info << "Creating field flameSens" << nl << endl;
    volScalarField flameSens
    (
        IOobject
        (
            "flameSens",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("flameSens", dimless, 0.0)
    );
    
    //- Thickened flame properties are read in from a file 
    IOdictionary thickenedFlameProperties
    (
        IOobject
        (
            "thickenedFlameProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //- Choose constant or dynamic thickening
    Switch dynamicTh = thickenedFlameProperties.lookupOrDefault<Switch>("dynamicTh", true);

    //- Choose local or global thickening
    Switch localTh = thickenedFlameProperties.lookupOrDefault<Switch>("localTh", true);
    
    //- Create LESdelta to compute thickening factor
    Info << "Creating field LESdelta" << nl << endl;
    volScalarField LESdelta
    (
        IOobject
        (
            "LESdelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("LESdelta", dimLength, 0.0005)
    );

    forAll(LESdelta, cellI)  // Loop over every cell of the mesh
    {  
         LESdelta[cellI] = Foam::pow(mesh.V()[cellI], 1.0/3.0);
    }	
    LESdelta.correctBoundaryConditions();
    
    //- Create maximum thickening factor in every cell
    volScalarField TFmax
    (
       IOobject
       (
           "TFmax",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionedScalar("TFmax", dimless, 1.0)
    );
    if (localTh) 
    {
        scalar DeltaEmax = readScalar(thickenedFlameProperties.lookup("DeltaEmax"));
        dimensionedScalar globalMaxDelta("globalMaxDelta", dimLength, DeltaEmax);
        TFmax = max(1.0,LESdelta/globalMaxDelta);
    } 
    else 
    {
        scalar TFmaxScalar = readScalar(thickenedFlameProperties.lookup("maxthickeningFactor"));
        TFmax = TFmax * TFmaxScalar;
    }

    //- Create thickening factor in every cell
    Info << "Creating field TF" << nl << endl;
    volScalarField TF
    (
        IOobject
        (
            "TF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TF", dimless, 1.0)
    );
    
//- Create function to compute efficiency function

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

    //- Create uPrime for computing efficiency function
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

    //- Create laminar flame velocity
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

    //- Create laminare flame thickness
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

    //- Create equivalence ratio
    Info << "Creating field eqRatio" << endl;
    volScalarField eqRatio
    (
        IOobject
        (
            "eqRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("eqRatio", dimless, 1.0)
    );
    
    scalar eqRatioConst = readScalar(thickenedFlameProperties.lookup("eqRatioConst"));
    eqRatio = eqRatio * eqRatioConst;
    Info << "Equivalence ratio = " << eqRatioConst << nl << endl;
    
    //- Create efficiency function
    autoPtr<efficiencyFunction> efficiency
    (
        efficiencyFunction::New(U, TF, eqRatio, gammafit, uPrime, dynsl0, dyndeltal0)
    );
    
  //  scalar effConst(thickenedFlameProperties.lookupOrDefault<scalar>("effConst", 1.0));
  //  Info << "Efficiency is multiplied by a constant = " << effConst << nl << endl;

//- Read properties for dynamic efficiency constant variation
    
    IOdictionary effFuncProperties
    (
        IOobject
        (
            "efficiencyFunctionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    scalar effConst(effFuncProperties.lookupOrDefault<scalar>("effConst", 1.0));
    
    Switch dynEffConst(effFuncProperties.lookupOrDefault<Switch>("dynEffConst", false));
    
    scalar effConstChangeRate(effFuncProperties.lookupOrDefault<scalar>("effConstChangeRate", 0.1));
    
    scalar maxPartProgVar(effFuncProperties.lookupOrDefault<scalar>("maxPartProgVar", 0.55));
    
    scalar crMaxPartProgVar(effFuncProperties.lookupOrDefault<scalar>("crMaxPartProgVar", 0.65));
    
    scalar minPartProgVar(effFuncProperties.lookupOrDefault<scalar>("minPartProgVar", 0.1));
    
    label minNumOfSteps(effFuncProperties.lookupOrDefault<label>("minNumOfSteps", 1000));

    
//- Create fields to evaluate the Laplacian on a coarser mesh
   
    //- Create field for velocity interpolated to the coarser mesh
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
   
    //- Create field for Laplacian operator on the coarser mesh
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
        dimensionedVector("coarseLap", dimensionSet(0,-1,-1,0,0,0,0), vector::zero)
    );
    
    //- Create field for Laplacian operator interpolated to the basic mesh 
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
    
 // Fields for flamelet closure
 
    //- Create the progress variable source term
/*    Info << "Creating field pvSource" << nl << endl;
    volScalarField pvSource
    (
        IOobject
        (
            "pvSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pvSource", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );

    //- Create the thickened progress variable source term
    Info << "Creating field pvSourceATF" << nl << endl;
    volScalarField pvSourceATF
    (
        IOobject
        (
            "pvSourceATF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pvSourceATF", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );*/
    
 // Fields and parameters for the turbulence forcing
 
    const dictionary turbForcing_(transportProperties_.subDict("turbForcing"));
    
    // Turbulence forcing control switch
    // true: turbulence forcing is on
    // false: turbulence forcing is off
    Switch forcedTurb(turbForcing_.lookupOrDefault<Switch>("forcedTurb", false)); 
    
    // Turbulent kinetic energy at the inlet
    dimensionedScalar k0("k0", dimensionSet(0,2,-2,0,0,0,0), readScalar(turbForcing_.lookup("k0")));
    
    // Turbulent velocity fluctuation at the inlet
    dimensionedScalar uprime0("uprime0", dimVelocity, readScalar(turbForcing_.lookup("uprime0")));
    
    // Turbulent length scale at the inlet
    dimensionedScalar lt0("lt0", dimLength, readScalar(turbForcing_.lookup("lt0")));
    
    // Number of domain cells in direction of flame propagation
    label Nx(readLabel(turbForcing_.lookup("Nx")));
    
    // Laminar flame speed
    scalar sL(readScalar(turbForcing_.lookup("sL")));
    
    // Desired flame position from syntheticTurbulenceInletUMeanControl
    scalar Xtarg(readScalar(turbForcing_.lookup("Xtarg")));
    
    // Maximal temperature in a respective laminar flame
    scalar TmaxL(readScalar(turbForcing_.lookup("TmaxL")));
    
    if (forcedTurb)
    {
        Info << "Turbulence forcing is switched on when t > " << Xtarg/sL << " s"<< nl
             << "Forced turbulence parameters:" << nl
             << "\t k0:                " << k0 << nl
             << "\t uprime0:           " << uprime0 << nl
             << "\t lt0:               " << lt0 << nl
             << "\t Nx:                " << Nx << nl
             << "\t TmaxL:             " << TmaxL << nl<< endl;       
    }
    else
    {
        Info << "Turbulence forcing is switched off" << nl << endl;
    }
    
    //- Create Force field for turbulent forcing
    volVectorField Force
    (
        IOobject
        (
            "Force",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Force", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    );
    
//////////////////////
// Chemistry fields //
//////////////////////

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

DimensionedField<scalar, volMesh> RRi
(
    IOobject
    (
        "RRi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("RRi", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);
