tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(rho*D, Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
              + fuelCloud.SYFuDNS(mesh,filter,Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);
            
            if (!pimple.finalIter())
            {

                YiEqn.solve(mesh.solver("Yi"));
            }
            else
            {
                YiEqn.solve(mesh.solver("YiFinal"));

            }

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
}
