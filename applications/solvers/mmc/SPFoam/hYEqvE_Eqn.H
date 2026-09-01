/*---------------------------------------------------------------------------*\
|                                hYEqvE Eqn                                   |
|-----------------------------------------------------------------------------|
| Solves species transport and species enthalpy equation                      |
\*---------------------------------------------------------------------------*/

// The coupling model lives on the particle cloud. With the PDF method switched
// off (or before any cloud has been created) there is nothing to ask, so fall
// back to plain transport: every species is solved and the relaxation source is
// switched off through a zero Indicator. XiEqn.H guards pSets the same way.
const bool pdfCoupled = pdfMethodOn && pSets.size();

const boolList solveEqvSpecie
(
    pdfCoupled
  ? pSets.first().coupling().solveEqvSpecie()
  : boolList(YEqvE.size(), true)
);

forAll(YEqvE, i)
{
    if (solveEqvSpecie[i])
        fieldsY.add(YEqvE[i]);
}

fieldsY.add(hEqvE);

tmp<fv::convectionScheme<scalar> > mvConvectionYi
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fieldsY,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

scalar dc0 = runTime.elapsedCpuTime();

if (pdfCoupled && pimple.firstIter())
{
    pSets.first().coupling().EqvETargetValues(Xi,YEqvETarget,TEqvETarget);

    //- This Step could be avoided to speed up
    if (pSets.first().coupling().enabled())
        hEqvETarget = thermoEqvETarget.he(p,TEqvETarget);

    TTARGET = TEqvETarget;
}

scalar dce = runTime.elapsedCpuTime() - dc0;

Info << "Execution Time for dens coupling= " << dce << " s" << endl;

label inertIndex = -1;

volScalarField Yt(0.0*YEqvE[0]);

// Zero stand-in used when there is no cloud to provide a coupling indicator
autoPtr<volScalarField> noIndicatorPtr;

if (!pdfCoupled)
{
    noIndicatorPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "IndicatorNoPDF",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
}

volScalarField& Indicator =
    pdfCoupled ? pSets.first().coupling().Indicator() : noIndicatorPtr();

forAll(YEqvE, specieI)
{
    if (YEqvE[specieI].name() != inertSpecie)
    {
        volScalarField& Yi = YEqvE[specieI];

        volScalarField& YTarget = YEqvETarget[specieI];

        if (!solveEqvSpecie[specieI])
        {
            Yi.primitiveFieldRef() = 0.0;
        }
        else
        {

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
//              + fvm::div(phi, Yi, "div(phi,Yi_h)")
              + mvConvectionYi->fvmDiv(phi, Yi)
              - fvm::laplacian(rho * DEff, Yi)
                ==
                (rho * (YTarget - Yi) / tauRel) * Indicator
            );

            YiEqn.relax();

            if (!pimple.finalIter())
            {

                YiEqn.solve(mesh.solver("Yi"));
            }
            else
            {
                YiEqn.solve(mesh.solver("YiFinal"));

            }

            Yi.max(0.0);

            Yi.min(1.0);

            Yt += Yi;
        }
    }
    else
    {
        inertIndex = specieI;
    }
}

if (inertIndex == -1)
{
    FatalErrorInFunction
        << "inertSpecie " << inertSpecie << " from the thermophysical "
        << "properties was not found among the transported species." << nl
        << "Without it the mass fractions cannot be closed."
        << exit(FatalError);
}

YEqvE[inertIndex] = scalar(1) - Yt;
YEqvE[inertIndex].max(0.0);

{
    fvScalarMatrix hEqvE_Eqn
    (
        fvm::ddt(rho, hEqvE)
//      + fvm::div(phi, hEqvE, "div(phi,Yi_h)")
      + mvConvectionYi->fvmDiv(phi, hEqvE)
      - fvm::laplacian(rho * DEff, hEqvE)
        ==
        (rho * (hEqvETarget - hEqvE) / tauRel) * Indicator
    );

    hEqvE_Eqn.relax();

    hEqvE_Eqn.solve();
}
// force the saving of the old-time values
// had been removed in OpenFOAM v2012, see issue #21
//psi.oldTime();
//thermo.correct();

volScalarField& T = thermo.T();

hEqvE.correctBoundaryConditions();
T.correctBoundaryConditions();
p.correctBoundaryConditions();
rho.correctBoundaryConditions();
thermo.correct(); 

rho = thermo.rho();
