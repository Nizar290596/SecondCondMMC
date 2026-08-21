Info << "==========================================" << endl;
Info << "Create filtered fields - START            " << endl;
Info << "==========================================" << endl;

// If the flag filterDNSFields is set to false, the original DNS velocity,
// density, temperature etc. is used. 
// Further, all fields are constructed and registered in the DNS mesh!
// To allow for a general structure refPtr of the fields are used and the 
// references are set accordingly.

refPtr<volScalarField> rhoFilteredRef;
refPtr<volVectorField> UFilteredRef;
refPtr<volVectorField> gradRhoFilteredRef;
refPtr<mmcVarSet> XiFilteredRef;

Info << "Create filtered turbulence model..." << endl;
// Generate the filtered turbulence properties
filteredTurbulence turbulenceFiltered
(
    mesh,
    filterMesh,
    rho,
    U,
    filterDNSFields,
    Sc_,
    Sct_
);

turbulenceFiltered.correct
(
    U,
    rho,
    p,
    T,
    Y
);

// Create filtering object
const filterDNSField& filter = turbulenceFiltered.filter();

// References for easier field handling
const volScalarField& DFiltered = turbulenceFiltered.D();
const volScalarField& DtFiltered = turbulenceFiltered.Dt();
const volScalarField& DEffFiltered = turbulenceFiltered.DEff();

if (filterDNSFields)
{
    rhoFilteredRef.cref(turbulenceFiltered.rho());
    UFilteredRef.cref(turbulenceFiltered.U());

    // Create the filtered mmcVarSet
    XiFilteredRef.reset
    (
        new mmcVarSet(entries,filterMesh)
    );

    Info << "Creating field gradRho" << nl << endl;
    gradRhoFilteredRef.reset
    (
        new volVectorField
        (
            IOobject
            (
                "gradRho",
                runTime.timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(rhoFilteredRef.cref() * DEffFiltered)
        )
    );
}
else
{
    rhoFilteredRef.cref(rho);
    UFilteredRef.cref(U);
    gradRhoFilteredRef.reset
    (
        new volVectorField
        (
            IOobject
            (
                "gradRho",
                runTime.timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(rho * DEffFiltered)
        )   
    );
    XiFilteredRef.cref(Xi);
}

// For easier handling provide const references here
const volScalarField& rhoFiltered = rhoFilteredRef.cref();
const volVectorField& UFiltered = UFilteredRef.cref();
const volVectorField& gradRhoFiltered = gradRhoFilteredRef.cref();
const mmcVarSet& XiFiltered = XiFilteredRef.cref();



