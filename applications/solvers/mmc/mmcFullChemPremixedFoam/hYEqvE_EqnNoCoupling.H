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

reaction->correct();
Qdot = reaction->Qdot();

label inertIndex = -1;

volScalarField Yt(0.0*YEqvE[0]);

forAll(YEqvE, specieI)
{
    if (YEqvE[specieI].name() != inertSpecie)
    {
        volScalarField& Yi = YEqvE[specieI];

            auto reactionMatrix = reaction->R(Yi);
            RRi = reactionMatrix & Yi;

            fvScalarMatrix YiEqn // CHANGE THIS
            (
                fvm::ddt(rho, Yi)
              + mvConvectionYi->fvmDiv(phi, Yi)
              - fvm::laplacian(TF * effConst * efficiency->eff() * rho * D + (1.0 - flameSens) * rho * Dt, Yi)
                ==
                effConst * RRi  * efficiency->eff()/ TF
                + fvOptions(rho, Yi)
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
        
            Yi.max(0.0);
            
            Yi.min(1.0);
        
            Yt += Yi;
//        }
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
      + mvConvectionYi->fvmDiv(phi, hEqvE)
      - fvm::laplacian(TF * effConst * efficiency->eff() * rho * D + (1.0 - flameSens) * rho * Dt, hEqvE)
        ==
        effConst * Qdot * efficiency->eff()/ TF
        + fvOptions(rho, hEqvE)
    );

    hEqvE_Eqn.relax();
    
    fvOptions.constrain(hEqvE_Eqn);
    
    hEqvE_Eqn.solve();
       
    fvOptions.correct(hEqvE);
}

psi.oldTime();
thermo.correct();

Info<< "min/max(T) = "
    << min(T).value() << ", " << max(T).value() << endl;
    
Info<< "min/max(Qdot) = "
    << min(Qdot).value() << ", " << max(Qdot).value() << endl;

