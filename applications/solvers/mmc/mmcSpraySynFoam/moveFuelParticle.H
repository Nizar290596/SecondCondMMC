if(fuelParticleOn)
{
    Info << "Moving Fuel Particles from solver" << endl;
 
    forAll(fSets,pSi)
    {
        fSets[pSi].evolve();
	    dStats.particleSampling(fSets[pSi]);
    }    
}
