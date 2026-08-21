if(true){
    solve(fvm::ddt(inflowU)  - fvm::laplacian(nu,inflowU)); // Diffusion, constant nu
    divU=fvc::div(inflowU);
    var.normalise();                                // re-normalise variances
    var.corTensor();
}
else{
    Info << "Using PISO Loop" << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(inflowU)
      + fvm::div(phi, inflowU)
      - fvm::laplacian(nu,inflowU)
    );

    if (true)                            // momentum predictor
    {
        solve(UEqn == -fvc::grad(inflowP));
    }

    // --- PISO loop
    for (int corr=0; corr<2; corr++)
    {
        volScalarField rAU(1.0/UEqn.A());
        
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        volVectorField HbyA("HbyA", inflowU); 
               
        HbyA = rAU*UEqn.H();
        
        const volScalarField one(UEqn.A()/UEqn.A());

        phi = (fvc::interpolate(HbyA) & auxMesh.Sf())
            + rAUf*fvc::ddtCorr(one, inflowU, phi);

        adjustPhi(phi, inflowU, inflowP);

        // Pressure corrector
        fvScalarMatrix pEqn
        (
           fvm::laplacian(rAU, inflowP) == fvc::div(phi)
        );

        pEqn.setReference(0, 0);
        pEqn.solve();
        phi -= pEqn.flux();

      //  #include "continuityErrs.H"
        inflowU -= rAU*fvc::grad(inflowP);
        inflowU.correctBoundaryConditions();
    }
    var.normalise2();                                // re-normalise variances
    var.corTensor();
}
