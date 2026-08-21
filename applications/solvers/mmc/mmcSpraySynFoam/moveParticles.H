/*---------------------------------------------------------------------------*\
|                               Move Particles                                |
|-----------------------------------------------------------------------------|
| Move particles and solve the reaction.                                      |
\*---------------------------------------------------------------------------*/

if(pdfMethodOn)
{
    Info << "Moving Particles" << endl;
 
    forAll(pSets,pSi)
    {
        Info << "\tParticle mixing..." << endl;
        
        basicAerosolReactingPopeParticle::trackingData td(pSets[pSi]);

        // Particle injection at inflow boundaries
        pSets[pSi].inflowBoundary().inflow();

        pSets[pSi].synthesis().synthesize();

        // Report time required for mixing
        const scalar t0mix = runTime.elapsedCpuTime();

        pSets[pSi].mixing().Smix();

        const scalar t1mix = runTime.elapsedCpuTime();

        Info << "\tFinish mixing in " << t1mix - t0mix<< "s" <<endl;



        Info << "\tSolve reaction..." << endl;
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
        Info << "\tFinish solve reaction in " << t1solve - t0solve<< "s" <<endl;
        
        pSets[pSi].composition().correctMassFractions();
        
        // Sample particle statistics for this cloud 
        pStats.particleSampling(pSets[pSi]);
        
        Info << "DONE"<<endl;
    }
}
