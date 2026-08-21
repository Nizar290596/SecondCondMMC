#include "createRDeltaT.H"


// =============================================================================
// Create Thermo Model for Eulerian DNS Field
// =============================================================================

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");



basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

// We need the density field for the solution of the mass conservation equation
// rhoEqn.H. The result of rhoEqn will be used for the momentum transport!
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
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

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;
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

// =============================================================================
// Construct Fields for mmc cloud
// =============================================================================
// Note: Here all fields are first created on the Eulerian DNS mesh. Only in 
//       the subsequent step of createFilteredFields.H the filtered fields
//       on the LES mesh are constructed. 


// Print the mmcFoam welcome message and construct the mmcVarSet
// =============================================================

//- Define location of file with mmcVariables definitions
const fileName solverDir = mmcStatusMessage::solverBasePath("applications/solvers/mmcDNS/mmcDNSFoam"); 

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


// We need the transport properties dictionary to read the turbulent and laminar
// Schmidt number. Standard compressible OpenFOAM solvers do not have a
// transportProperties file, but also do not know a Schmidt number...
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

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsXi;

forAll(Xi, i)
{
    fieldsXi.add(Xi.Vars(i).field());
}

// Laminar diffusivity on the DNS mesh
volScalarField D
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()/Sc_
);


// =============================================================================
// Construct Fields for mmc cloud -- END
// =============================================================================

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());


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

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"

