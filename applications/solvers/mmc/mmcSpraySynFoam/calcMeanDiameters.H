Info<<"Calculating D10(AMD) and D32(SMD) variables for FP"<<endl;

    // Initialse the fields at every time step
    forAll(numAMD,i)
    {
        numAMD[i] = 0.0;   // D10 numerator
        denAMD[i] = 0.0;   // D10 denominator
        
        numSMD[i] = 0.0;   // D32 numerator
        denSMD[i] = 0.0;   // D32 denominator
    }

    //- Accumulate the varaibles in each computational cell
    forAllIters(fSets[0], iter)
    {	
        const label celli = iter().cell();

        // Calculating D10 numerators and denominators
        label p = 1;
        label q = 0;
        numAMD[celli] += iter().nParticle()*pow(iter().d(), p);
        denAMD[celli] += iter().nParticle()*pow(iter().d(), q);

        // Calculating D32 numerators and denominators
        label m = 3;
        label n = 2;
        numSMD[celli] += iter().nParticle()*pow(iter().d(), m);
        denSMD[celli] += iter().nParticle()*pow(iter().d(), n);
    }
