{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(rho*D, he) //in MMCLES, we use the rho * DEff,
     ==
        Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    // force the saving of the old-time values
    // had been removed in OpenFOAM v2012, see issue #21
    psi.oldTime(); // from MMCLES

    thermo.correct();
    
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

    Info<< "min/max(Qdot) = "
        << min(Qdot).value() << ", " << max(Qdot).value() << endl;
}
