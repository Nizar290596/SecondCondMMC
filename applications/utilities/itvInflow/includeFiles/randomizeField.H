// Randomization
if(args.found("nr")){Info << "skipping randomization" << endl;}    // dont randomize
else{
    Info << "Randomize: " ;
    Random ranGen_(time(nullptr));                  // initialize random generator
    scalar msize = auxMesh.V().size();
    for(int i=0; i<msize; i++){                     // add randomness and normalize
        vector input = ranGen_.globalSample01<vector>();          // it with cell volume Kempf "Efficient Generation ..." p.74
        dimensionedScalar V(auxMesh.V()[i]);
        for(int j=0;j<3;j++){
            inflowU[i][j] = (input[j]-0.5)/(sqrt(V).value());// normalize with cell volume
        }
    }
    Info << " done" << endl;
}
