/*---------------------------------------------------------------------------*\
|                             Read Flamelet Table                             |
|-----------------------------------------------------------------------------|
|              Reads the flamelet table and its properties                    |
| One dimentional version (for 2-dimentional tables check mmcPremixedAtfFoam) |
\*---------------------------------------------------------------------------*/
#include "RectangularMatrix.H"

//- Reading properties associated with the flamelet table
Info << "Reading flameletTableProperties" << nl << endl;
IOdictionary flameletTableProperties
(
    IOobject
    (
        "flameletTableProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

fileName   FLTfile(fileName(flameletTableProperties.lookup("FLTFileName")).expand()); // Path to the famelet table file
label      nCol    = readScalar(flameletTableProperties.lookup("nCol")); // no. the famelet table columns (1...nCol)
label      nRow    = readScalar(flameletTableProperties.lookup("nRow")); // no. the famelet table rows (1...nRow)
label      pvDisc  = readScalar(flameletTableProperties.lookup("pvDisc"));  // no. of discretization points for progress variable (needed if tables contain progress var and mixture fraction)
label      sourceI = readScalar(flameletTableProperties.lookup("pvSourceI")); // index of a column containing the progress variable source term
label      YCO2I = readScalar(flameletTableProperties.lookup("YCO2I")); // index of a column containing the CO2 mass fraction

//- Reading the flamelet table
Info << "Reading flamelet table" << nl << endl;

using std::ifstream;
ifstream FlameletTable;
FlameletTable.open(FLTfile.c_str());

RectangularMatrix <scalar> FGMtable(nRow,nCol);

for (label i = 0; i < nRow; i++)
{
    for (label j = 0; j < nCol; j++) 
    {
        FlameletTable >> FGMtable[i][j];
    }
}

FlameletTable.close();

//- Get the equilibrium value of the CO2 mass fraction for conditioning variable normalizaation
scalar YCO2eq = 1e30;
YCO2eq = FGMtable[pvDisc][YCO2I];
Info << "Equilibrium value of the CO2 mass fraction, YCO2eq = "<< YCO2eq << nl << endl;
