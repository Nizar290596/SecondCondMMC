/*---------------------------------------------------------------------------*\
|          Calculate turbulence forcing and flow velocity                     |
|-----------------------------------------------------------------------------|
| Turbulence forcing is calculated if velocity fluctuations from the inlet    |
| have reached the flame position                                             |
| UEqn is solved to get the velocity field                                    |
\*---------------------------------------------------------------------------*/

// Calculate turbulence forcing if fluctuations from the inlet reached the flame position (runTime.value() > Xtarg/sL)
if (forcedTurb == true && runTime.value() > Xtarg/sL)
{
    Info << "Solving Turbulence Forcing" << endl;
    vector rhoUx[Nx];
    vector Ux[Nx];
    scalar kx[Nx];
    scalar rhox[Nx];            
    scalar count[Nx];

    for (label j = 0; j < Nx; j++)
    {
        rhoUx[j] = vector::zero;
	    Ux[j] = vector::zero;
	    kx[j] = 0;
	    rhox[j] = 0;
        count[j] = 0;
    }
 
    forAll(mesh.cells(), cellI)
    {
	    scalar X =  mesh.cellCentres()[cellI].component(0);
        scalar dx = Foam::pow(mesh.V()[cellI], 1.0/3.0);
	    
	    for (label i = 0; i < Nx; i++)
        {
	        if (X < (i+1)*dx && X >= i*dx)
	    	{
			    rhoUx[i] = rhoUx[i] + U[cellI]*rho[cellI];
			    rhox[i] = rhox[i] + rho[cellI];
			    count[i] = count[i] + 1;
            }
        }
    }
 
    for (label i = 0; i < Nx; i++)
    {
	    reduce(rhoUx[i], sumOp<vector>());
	    reduce(rhox[i], sumOp<scalar>());
	    reduce(count[i], sumOp<scalar>());
    }
 
    for (label i = 0; i < Nx; i++)
    {
	    rhoUx[i] = rhoUx[i] / count[i];
	    rhox[i] = rhox[i] / count[i];
	    Ux[i] = rhoUx[i] / rhox[i];
    }

    forAll(mesh.cells(), cellI)
    {
	    scalar X =  mesh.cellCentres()[cellI].component(0);
        scalar dx = Foam::pow(mesh.V()[cellI], 1.0/3.0);
        
	    for (label i = 0; i < Nx; i++)
	    {
	    	if (X < (i+1)*dx && X >= i*dx)
	    	{
			    kx[i] = kx[i] + 1.0/2.0 * magSqr(U[cellI] - Ux[i]) * rho[cellI];
            }
	    }
    }
 
    for (label i = 0; i < Nx; i++)
    {
	    reduce(kx[i], sumOp<scalar>());
    }
 
    for (label i = 0; i < Nx; i++)
    {
	    kx[i] = kx[i] / count[i] / rhox[i];
    }

    dimensionedScalar A = uprime0/lt0/3.0;

    forAll(mesh.cells(), cellI)
    {
        scalar X =  mesh.cellCentres()[cellI].component(0);
        scalar dx = Foam::pow(mesh.V()[cellI], 1.0/3.0);
        
	    for (label i = 0; i < Nx; i++)
	    {
	        if (X < (i+1)*dx && X >= i*dx)
	    	{
		        if (thermo.T()[cellI] < thermo.T().boundaryField()[mesh.boundaryMesh().findPatchID("INLET")][0] + 0.01*(TmaxL - thermo.T().boundaryField()[mesh.boundaryMesh().findPatchID("INLET")][0]))
		        {
			        Force[cellI] = A.value() * k0.value() / kx[i] * ( rho[cellI] * U[cellI] - rhoUx[i] );
		        }
            }
	    }
    }
} 
Force.correctBoundaryConditions();

// Solve the Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g     
      + Force    // Forcing is added
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        K = 0.5*magSqr(U);
    }
