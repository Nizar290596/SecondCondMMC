 

    Info<<"Calculating Volume flux for liquid fuel particle"<<endl;

    //- set all flux values zero at begining of time step
    forAll(particleFlux,i)
    {
        evaporationRate[i] = 0.0;
        particleFlux[i]    = 0.0;
        particleFlux10[i]  = 0.0;
        particleFlux20[i]  = 0.0;
        particleFlux30[i]  = 0.0;
        particleFlux40[i]  = 0.0;
        particleFlux50[i]  = 0.0;
    }	        


    //- Accumulate all parcels mass in a cell and divide by rho
    forAllIters(fSets[0], iter)
    {
		
		const label cellI = iter().cell();			
		particleFlux[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
        evaporationRate[cellI]+=iter().mFlux() * iter().d() * iter().d() * constant::mathematical::pi;

        if(iter().d()<10e-06)
        {
            particleFlux10[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
        }
        else
        {
            if(iter().d()<20e-06)
            {
                particleFlux20[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
            }
            else
            {
                if(iter().d()<30e-06)
                {
                    particleFlux30[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
                }
                else
                {
                    if(iter().d()<40e-06)
                    {
                        particleFlux40[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
                    }
                    else
                    {
                        if(iter().d()<50e-06)
                        {
                            particleFlux50[cellI]+=iter().nParticle()*iter().mass()*iter().U().z()/iter().rho();
                        }
                    }
                }
            }
        }
        
    }

    //- Calculate volume flux 
    forAll(mesh.V(),i)
    {
        particleFlux[i]   /= mesh.V()[i];
        particleFlux10[i] /= mesh.V()[i];
        particleFlux20[i] /= mesh.V()[i];
        particleFlux30[i] /= mesh.V()[i];
        particleFlux40[i] /= mesh.V()[i];
        particleFlux50[i] /= mesh.V()[i];
    }

	

