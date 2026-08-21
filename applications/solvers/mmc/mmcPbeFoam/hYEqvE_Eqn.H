
forAll(YEqvE, i)
{
    if (pSets.first().coupling().solveEqvSpecie()[i])
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

if (pimple.firstIter())
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

volScalarField& Indicator = pSets.first().coupling().Indicator();

forAll(YEqvE, specieI)
{
    if (YEqvE[specieI].name() != inertSpecie)
    {
        volScalarField& Yi = YEqvE[specieI];

        volScalarField& YTarget = YEqvETarget[specieI];

        if (!pSets.first().coupling().solveEqvSpecie()[specieI])
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

thermo.correct();

