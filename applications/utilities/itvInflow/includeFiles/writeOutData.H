Info << nl << "writing out planes" << endl;                // TODO replace by propper class
initializeOutput(runTime,inflowProp.patchName(),pp);

meshSearch meshsearch(auxMesh);

 const label skipCellsWO(max(1,inflowProp.intL()/(4*statisticProps.deltaZ())));    // to reduce amount of output data only 4 Values per structure are written

// debuging by zhijie, output all the slices of the whole domain 
//const label skipCellsWO(1);


scalar lambda(inflowProp.lambda());                                                        // scaling factor
const scalar deltaT(statisticProps.deltaZ()/lambda/inflowProp.ubulk()*skipCellsWO);// time step for write out. distance of layers /  bulk velocity
runTime.setDeltaT(deltaT);
runTime.setTime(runTime.startTime().value()-deltaT+timeOffset,0);        // start at startTime-deltaT to write out at startTime !! fishy !!
Info <<  "deltaT: " << deltaT << " seconds \n"<< endl;
label Oiter(0);

point position;
scalar Zpos;

autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New    // construct UInterp
(
        "cell",
        inflowU
);

while (runTime.loop() && !restart)
{
    Info<< "\r Time = " << runTime.timeName();
    fvPatchField<vector> boundaryBackup(origUField);

     for(int Iiter=0; Iiter<origUField.patch().Cf().size(); Iiter++){                 // add fluctuations to original patch
         Zpos = Oiter*skipCellsWO*statisticProps.deltaZ()+statisticProps.deltaZ()/2;
         position = origUField.patch().Cf()[Iiter]*lambda;                            //
         position[2] = 0;                                                            // set z pos to zero
         position += vector(0.5,0.5,Zpos);
         label cell     = meshsearch.findCell(position);                                // dont use the expensive mesh search
         origUField[Iiter] += UInterp().interpolate(position,cell);
     }

     writeOutPatch(runTime, inflowProp.patchName(), origUField);    // write out data
     origUField == boundaryBackup;                                // restore patch
     Oiter++;

    if( Oiter > statisticProps.nCellsZ()/skipCellsWO-1)
    {
        Info << "\nrestarting inflow data generation after: "  <<  runTime.timeName()  << " seconds \n";
        Oiter = 0;
        restart = true;
        timeOffset=runTime.timeOutputValue();
    }
}
Info << "\n\ndone" << endl;
