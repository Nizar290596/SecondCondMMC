void initializeOutput(const Foam::Time &runTime,const word clonePatchName, const polyPatch& pp){
    fileName sampleBD(runTime.constant()/"boundaryData");                    // dont hard code the patch name
    fileName samplePath(sampleBD/clonePatchName);
    mkDir(samplePath);
    fileName pointsFile(samplePath/"points");
    OFstream points(pointsFile);

    //points << "FoamFile\n{\nversion\t2.0;\nformat\tascii;\nclass\tvectorField;\nobject\tpoints;\n}\n"; // write stupid header
    //do not write the headers. to match the 5.x format.

    points << pp.faceCentres();                                                // write out face center points

}


inline void writeOutPatch(const Time &runTime,const word clonePatchName, const fvPatchField<vector> &U){
    // write out data
    fileName sampleBD(runTime.constant()/"boundaryData");
    fileName samplePath(sampleBD/clonePatchName);
    fileName timeDir(samplePath/runTime.timeName());
    mkDir(timeDir);
    fileName sampleFile(timeDir/"U");
    OFstream samples(sampleFile);
    //samples << "FoamFile\n{\nversion\t2.0;\nformat\tascii;\nclass\tvectorAverageField;\nobject\tvalues;\n}\n\n(0 0 0)\n"; // write stupid header
    //do not write the headers. to match the 5.x format.
    samples << "\n"<<  U.size() << "\n(" << endl;
    forAll(U,iter){
        samples << U[iter] << "\n";
    }
    samples << "\n)" << endl;
}
