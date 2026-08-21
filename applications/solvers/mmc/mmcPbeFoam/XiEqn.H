dimensionedScalar tauRel("tauRel", dimTime, 0);

label deltaIndices = mesh.time().timeIndex() - mesh.time().startTimeIndex() - 1;

tauRel.value() = max(pSets.first().coupling().tauRelaxStart().value() - deltaIndices*pSets.first().coupling().tauRelaxDelta(), pSets.first().coupling().tauRelaxTarget().value());

if (pSets.first().coupling().tauUnits() == "timestep")
{
    tauRel.value() = tauRel.value()*mesh.time().deltaT().value();
}

if (pSets.first().coupling().tauRelaxBlending())
{
    Info << "tauRealaxStart = " << pSets.first().coupling().tauRelaxStart().value() << endl;
    Info << "tauRelaxTarget = " << pSets.first().coupling().tauRelaxTarget().value() << endl;
    Info << "tauRelaxDelta = " << pSets.first().coupling().tauRelaxDelta() << endl;
    Info << "tauRealax = " << tauRel.value() << endl;
}

if (tauRel.value() < mesh.time().deltaT().value())
    tauRel.value() = mesh.time().deltaT().value();
    
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
    volScalarField& XiiTarget = XiTarget.Vars(XiI).field();
    
    if (isA<couplingVar>(Xi.Vars(XiI)) || (Xi.Vars(XiI).refType()=="interpolated"))
    {
        
        if (Xi.Vars(XiI).passive())
            XiiTarget = Xii;            
        else
            XiiTarget = 0.0*Xii; // Need to add a way to calculate targets or 
                                 // a source term
            
        fvScalarMatrix XiiEqn
        (
            fvm::ddt(rho,Xii)
//          + fvm::div(phi,Xii,"div(phi,Xii)")
          + mvConvectionXi->fvmDiv(phi, Xii)
          - fvm::laplacian(rho * DEff,Xii)
            ==
            rho * (XiiTarget - Xii) / tauRel
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
