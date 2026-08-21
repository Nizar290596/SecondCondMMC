/*---------------------------------------------------------------------------*\
|                      Read Mixing Time Scale Table                           |
|-----------------------------------------------------------------------------|
| Reads tabulated mixing time scale for laminar premixed flames (tau_LPF)     |
| as a function of the progress variable                                      |
\*---------------------------------------------------------------------------*/

Info << "Reading tabulated LPF mixing time scale:" << nl;
 
std::ifstream MixTimeTable;
MixTimeTable.open("mixTimeTab.dat");
    
label nRowsMixTime = 501; 
RectangularMatrix <scalar> mixingTime(nRowsMixTime, 2);
scalar mixTimeLPFmin = 1e30;

for (label i = 0; i < nRowsMixTime; i++)
{
    for (label j = 0; j < 2; j++)
    {
        MixTimeTable >> mixingTime[i][j];
    }
    
    Info << "c: " << mixingTime[i][0] << token::TAB << " mixTime LPF: " << mixingTime[i][1] << nl;
    
    if (mixTimeLPFmin > mixingTime[i][1]) mixTimeLPFmin = mixingTime[i][1]; // Get minimal LPF mixing time scale value to calculate the substeps for particle mixing
}

MixTimeTable.close();
