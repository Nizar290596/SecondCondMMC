/*---------------------------------------------------------------------------*\
|                               Move Particles                                |
|-----------------------------------------------------------------------------|
| Move particles and solve the reaction.                                      |
| Sample particle statistics if activated in the cloud dictionary             |
|                                                                             |
| Timestep sequence                                                           |
|   1. inflowBoundary().inflow()          particle injection                  |
|   2. mixing().Smix()                    first conditioning (all particles)  |
|   3. updatePhiReaction(Δt)              W(φ) source  (all particles)        |  <- guarded
|   4. updateOUProcess(Δt)               OU advance + φ° update (subset)     |  <- guarded
|   5. secondCondMixing().Smix()          second conditioning (subset)        |  <- guarded
|   6. solve()                            Lagrangian transport + chemistry     |
|   7. Sreact()                           load-balanced chemistry              |
|                                                                             |
| Steps 3–5 are only executed when secondCondMixingEnabled() is true.        |
\*---------------------------------------------------------------------------*/

if(pdfMethodOn)
{
    Info << "Evolve Particles -- Mixing, Moving, and Reaction" << endl;
 
    forAll(pSets,pSi)
    {
        Info << "\tParticle mixing..." << endl;
        
        basicReactingPopeParticle::trackingData td(pSets[pSi]);

        // Particle injection at inflow boundaries
        pSets[pSi].inflowBoundary().inflow();

        // Report time required for mixing
        const scalar t0mix = runTime.elapsedCpuTime();

        pSets[pSi].mixing().Smix();

        const scalar t1mix = runTime.elapsedCpuTime();

        Info << "\tFinish mixing in " << t1mix - t0mix<< "s" <<endl;

	// Steps 3–5: second-conditioning block.
        // secondCondEnabled is true when the secondConditioning sub-dictionary
        // is present, enabled = true, and the model pointer was constructed.
        // updatePhiReaction (step 3) applies to ALL particles — φ is a
        // first-level quantity — but it only exists when second conditioning
        // is part of the setup, so it shares the same guard as steps 4 and 5.
        // When second conditioning is off the solver goes directly to step 6,
        // identical to the original baseline solver.
        const bool secondCondEnabled = pSets[pSi].secondCondMixingEnabled();
 
        if (secondCondEnabled)
        {
            const scalar deltaT = runTime.deltaTValue();
 
            // 3. W(φ) reaction — all particles, after Smix() has completed
            //    for the whole cloud so φ values are post-mixing
            pSets[pSi].updatePhiReaction(deltaT);
 
            // 4. Advance ω_OU and recompute φ° — flagged particles only,
            //    after reaction so φ° reflects the reacted φ
            pSets[pSi].updateOUProcess(deltaT);
 
        Info << "\tFinish mixing in " << t1mix - t0mix<< "s" <<endl;
            // 5. Second-conditioning mixing — uses updated φ° values
            //    (stored as eulerianFieldData.XiR()[0] in buildParticleList)
            const scalar t0sc = runTime.elapsedCpuTime();
 
            pSets[pSi].secondCondMixing().Smix();
 
            const scalar t1sc = runTime.elapsedCpuTime();
            Info << "\tFinish second-conditioning mixing in "
                 << t1sc - t0sc << "s" << endl;
        }

        Info << "\tMove particle and solve reaction..." << endl;
        const scalar t0solve = runTime.elapsedCpuTime();

	Pout<< "\t [solve] entering, local nParticles = " << pSets[pSi].size() << endl;

        // In solve() the particle move() function is called 
        // if balanceReactionLoad is true, solve(td) does not calculate chemistry
        pSets[pSi].solve(pSets[pSi],td);

	Pout<< "\t [solve] exited, local nParticles = " << pSets[pSi].size() << endl;

        // Sreact() is only defined in finiteRateChem, it does nothing for other reacting models
        if(pSets[pSi].balanceReactionLoad())
        {
            scalar t0react = runTime.elapsedCpuTime();

            pSets[pSi].reaction().resetCPUTimeAnalysis();//for analysis

	    Pout<< "\t [Sreact] entering" <<endl;

            pSets[pSi].reaction().Sreact();

	    Pout<< "\t [Sreact] exited" <<endl;

            pSets[pSi].reaction().writeCPUTimeAnalysis();//for analysis

            pSets[pSi].reaction().ISATCntIncrmnt();//timestep used by ISAT to clean tree etc.
                                                    // do nothing if finiteRate(ISAT) off
            
            scalar t1react = runTime.elapsedCpuTime();
            
            Info << "\tFinish balanced reaction in " << t1react- t0react<< "s" <<endl;
        }
        const scalar t1solve = runTime.elapsedCpuTime();
        Info << "\tFinish solve reaction and move in " << t1solve - t0solve<< "s" <<endl;
        
        pSets[pSi].composition().correctMassFractions();
        
        // Sample particle statistics for this cloud 
        pStats.particleSampling(pSets[pSi]);
        
        Info << "DONE"<<endl;
    }
    
}

