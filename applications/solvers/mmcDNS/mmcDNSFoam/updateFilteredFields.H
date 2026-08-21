Info << "Update filtered fields..." << endl;

// Always update the filteredTurbulence
turbulenceFiltered.correct
(
    U,
    rho,
    p,
    T,
    Y
);

if (filterDNSFields)
{
    // Update the reference fields
    forAll(Xi, j)
    {
        volScalarField& Xii = Xi.Vars(j).field();
        mmcVarSet& XiFiltered = XiFilteredRef.ref();
        volScalarField& XiiFiltered = XiFiltered.Vars(j).field();
        filter.updateField(Xii,rho,XiiFiltered);
    }
    gradRhoFilteredRef.ref() = fvc::grad(rhoFiltered*DEffFiltered);
}
else
{
    gradRhoFilteredRef.ref() = fvc::grad(rho*D);
}


