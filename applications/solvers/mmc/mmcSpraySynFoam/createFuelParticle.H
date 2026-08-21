    dictionary createFuelParticle = mesh.solutionDict().subDict("FUELPARTICLE");
    Switch fuelParticleOn = createFuelParticle.lookupOrDefault<Switch>("fuelParticleOn", false);

    label fSet = createFuelParticle.lookupOrDefault<label>("fSet",1);

    PtrList<basicSprayThermoCloud> fSets;
    fSets.clear();

    std::stringstream fName;
    
    volScalarField& f =Xi.Vars("f").field();

    // For the termo cloud of parcels OpenFOAM uses by default 
    // SLGThermo
    // Construct from main thermo model 
    SLGThermo slgThermo(mesh,thermo);

    if (!fuelParticleOn)
    {
        Info << "createFuelParticle method is switched off. \n" << endl;
    }
    else
    {
        Info << "createFuelParticle method is switched on Creating " << fSet << " set(s) of Fuel particles. \n" << endl;
        fName << "fuelCloud1";
        
        fSets.append
        (
            new basicSprayThermoCloud
            (
                fName.str(),
                rho,
                U,
                muEff,
                g,
                slgThermo,              
                f,                    	
                pSets.first(),
                magSqrGradf
            )
        );
    }
