/*---------------------------------------------------------------------------*\
|                      Stochastic Particles Evolution                         |
|-----------------------------------------------------------------------------|
| Move particles and solve the reaction, mix particles.                       |
| Sample particle statistics if activated in the cloud dictionary             |
\*---------------------------------------------------------------------------*/

if(pdfMethodOn)
{
    Info << "Evolve Particles -- Mixing, Moving, and Reaction" << endl;
    
    // Calculate detlaT for particle substeps (criterion: mixing extent < 0.5 for tau_LPF)
    scalar eulDeltaT = mesh.time().deltaT().value();
    
    scalar dtSubstep = eulDeltaT;
    
    label nSubsteps = 1;
    
    dictionary subModels(cloudProperties.subDict("subModels"));
    
    dictionary premixedMMCcurlCoeffs(subModels.subDict("premixedMMCcurlCoeffs"));
    
    //- Account (or not) for the ATF in the LPF mixing time scale distribution 
    //- true: tabulated LPF time scale is distributed over the entire reference space
    //- false: tabulated LPF is rescaled to act in a central section of the reference space with width equal to 1/TFmax
    const Switch ATFmixTimeOn(premixedMMCcurlCoeffs.lookup("ATFmixTime"));

    //- Calculate a number of time substeps for particles automatically or not 
    //- true: increase number of substeps until mixing extent < 0.5
    //- false: give number of substeps (numPartSubsteps) manually
    const Switch autoParticleTimeStep(premixedMMCcurlCoeffs.lookupOrDefault("autoParticleTimeStep", true));
    const label numPartSubsteps = premixedMMCcurlCoeffs.lookupOrDefault<label>("numPartSubsteps", 1);
    
    const scalar CL = readScalar(premixedMMCcurlCoeffs.lookup("CL"));
    
    
    if (autoParticleTimeStep)
    {
        //- Increase number of substeps until mixing extent < 0.5
        while(Foam::pow(dtSubstep, 0.5) >= (0.7 * CL * mixTimeLPFmin)) 
        {
            nSubsteps++;
            dtSubstep = eulDeltaT/nSubsteps;
        }
    }
    else 
    {
        nSubsteps = numPartSubsteps;
        dtSubstep = eulDeltaT/nSubsteps;
    }
 
    forAll(pSets,pSi)
    {        
        basicPremixedReactingPopeParticle::trackingData td(pSets[pSi]);

        // Particle injection at inflow boundaries
        pSets[pSi].inflowBoundary().inflow();
        
        Info << "\tTime step for the particles: "<< dtSubstep << " s, particle mixing and movement are divided in "<< nSubsteps << " substep(s)."<< endl;

        // Set deltaT to the substep value
        runTime.setDeltaT(dtSubstep, false);
        
        // Solve particle mixing and movement in substeps
        for (label i = 1; i <= nSubsteps; ++i) 
        {
            Info << "\tParticle evolution substep "<< i << endl;
            // Report time required for mixing
            const scalar t0mix = runTime.elapsedCpuTime();
            
            Info << "\tParticle mixing..." << endl;
            
            // Set tabulated LPF time scale for each particle based on the reference variable
            forAllIter(basicPremixedReactingPopeCloud, pSets[pSi], iter)
            {
                scalar cPrlxTab = 0.0;
                
                if (ATFmixTimeOn)
                {
                    scalar cPrlx = iter().XiR("cLES");
                    
                    interpolationCellPoint<scalar> TFmax_intp_(TFmax);
                    scalar TFmaxp = TFmax_intp_.interpolate(iter().position(),iter().cell()); 
                    
                    //- Rescaling of the reference space to actual particle flame thickness 
                    // assuming that TFmax value is approx. the same in the flame region (uniform mesh spacing)
                    if (cPrlx <= ((1 - 1/TFmaxp)/2))       { cPrlxTab = 0.0; }
                    else if (cPrlx >= ((1 + 1/TFmaxp)/2))  { cPrlxTab = 1.0; }
                    else
                    {
                         scalar TFp = 0.0;
                         
                         if (dynamicTh) { TFp = 1.0 + (TFmaxp - 1.0) * 16.0 * pow((cPrlx * (1.0 - cPrlx)), 2); }
                         else           { TFp = TFmaxp; }
                         
                         cPrlxTab = 0.5 + TFp * (cPrlx - 0.5);
                    }
                }
                else
                {
                    cPrlxTab = iter().XiR("cLES");
                }
                
                scalar pMixTime = 1e30;
                label rowI = 0;
            
                while(cPrlxTab >= mixingTime[rowI][0]) 
                {
                    if(cPrlxTab < mixingTime[rowI + 1][0]) 
                    {
                        pMixTime = mixingTime[rowI][1] + (mixingTime[rowI+1][1]-mixingTime[rowI][1]) * (cPrlxTab-mixingTime[rowI][0]) / (mixingTime[rowI+1][0]-mixingTime[rowI][0]);
                    }
                    rowI++;
                }
            
                iter().mixTime() = pMixTime * Foam::pow(dtSubstep, 0.5);
            }

            pSets[pSi].mixing().Smix();

            const scalar t1mix = runTime.elapsedCpuTime();

            Info << "\tFinish mixing in " << t1mix - t0mix<< "s" <<endl;

            Info << "\tMove particle and solve reaction..." << endl;
            const scalar t0solve = runTime.elapsedCpuTime();
            // if balanceReactionLoad is true, solve(td) does not calculate chemistry
            pSets[pSi].solve(pSets[pSi],td);

            // Sreact() is only defined in finiteRateChem, it does nothing for other reacting models
            if(pSets[pSi].balanceReactionLoad())
            {
                scalar t0react = runTime.elapsedCpuTime();

                pSets[pSi].reaction().resetCPUTimeAnalysis();//for analysis

                pSets[pSi].reaction().Sreact();

                pSets[pSi].reaction().writeCPUTimeAnalysis();//for analysis

                pSets[pSi].reaction().ISATCntIncrmnt();//timestep used by ISAT to clean tree etc.
                                                    // do nothing if finiteRate(ISAT) off
            
                scalar t1react = runTime.elapsedCpuTime();
            
                Info << "\tFinish balanced reaction in " << t1react- t0react<< "s" <<endl;
            }
            const scalar t1solve = runTime.elapsedCpuTime();
            Info << "\tFinish solve reaction and move in " << t1solve - t0solve<< "s" <<endl;
        
            pSets[pSi].composition().correctMassFractions();
        }
        
        // Set deltaT back to the initial (eulerian field) value
        runTime.setDeltaT(eulDeltaT, false);
        
        //- Calculate conditioning variable (progress variable on the particles) 
        Info << "\tUpdate reaction progress on particles " << endl;
               
        forAllIter(basicPremixedReactingPopeCloud, pSets[pSi], iter)
        {
            //- get index of CO2 mass fraction
            label idxCO2p = -1;
            forAll(pSets[pSi].composition().componentNames(), specieI)
            {
                if (pSets[pSi].composition().componentNames()[specieI] == "CO2")
                    idxCO2p = specieI;
            }
            //- normalize
            iter().XiC("cP") = iter().Y()[idxCO2p]/YCO2eq; 
        }
        
        // Sample particle statistics for this cloud 
        pStats.particleSampling(pSets[pSi]);
        
        Info << "DONE" << endl;
    }
    
}
