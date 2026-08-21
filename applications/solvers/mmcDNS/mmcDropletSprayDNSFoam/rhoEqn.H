/*---------------------------------------------------------------------------*\
|                               rhoEqn                                        |
|-----------------------------------------------------------------------------|
| Solve the density transport for consistent density                          |
\*---------------------------------------------------------------------------*/
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        fuelCloud.SrhoDNS(mesh,filter,rho)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);
}

