//  readLamTable.H: Read laminar properties in format
// "equivalenceRatio     laminarFlameThickness[mm]    laminarFlameSpeed[cm/s]" 

    using std::ifstream;

    Info << "Reading laminar properties table for efficiency calculation" << nl << endl;
 
    ifstream LamTable;

    LamTable.open(Lamfile.c_str());

    for (int i=0;i<nRow_;i++)
    {
        for (int j=0;j<nCol_;j++)
        {
            LamTable >> lamtable[i][j];
        }
    }

    LamTable.close();
