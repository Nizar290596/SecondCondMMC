
tmp<fv::convectionScheme<scalar> > mvConvectionXi
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fieldsXi,
        phi,
        mesh.divScheme("div(phi,Xii)")
    )
);

forAll(Xi, XiI)
{
    volScalarField& Xii = Xi.Vars(XiI).field();
    
    if (isA<couplingVar>(Xi.Vars(XiI)) || (Xi.Vars(XiI).refType()=="interpolated"))
    {
        fvScalarMatrix XiiEqn
        (
            fvm::ddt(rho,Xii)
          + mvConvectionXi->fvmDiv(phi, Xii)
          - fvm::laplacian(rho * D,Xii)
          ==
            fuelCloud.SrhoDNS(mesh,filter,Xii)
        );

        XiiEqn.relax();

        if (!pimple.finalIter())
        {   
            XiiEqn.solve(mesh.solver("Xii"));
        }
        else
        {
            XiiEqn.solve(mesh.solver("XiiFinal"));
        }

        Xii.max(0.0);
        
        Xii.min(1.0);
    }
    else
    {
        Info << "Internal Field is 0" <<endl;
        Xii.primitiveFieldRef() = 0.0;
    }
}

