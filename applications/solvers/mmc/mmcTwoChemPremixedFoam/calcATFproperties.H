/*---------------------------------------------------------------------------*\
|                             Calculate fields for ATF                        |
|-----------------------------------------------------------------------------|
| Calculates the thickening factor and the efficiency function                |
| Calculation of uPrime is done by first evaluating lapacian on a coarser mesh|
\*---------------------------------------------------------------------------*/

volScalarField& XiPV = Xi.Vars("cLES").field();

Info << "Calculating thickening factor for ATF" << endl;
forAll(XiPV,cellI)
{
    //- Calculate the flame sensor
    if (XiPV[cellI] >= 0.0 && XiPV[cellI] <= 1.0) 
    {
         flameSens[cellI] = 16.0 * pow((XiPV[cellI] * (1.0 - XiPV[cellI])), 2);
         if (flameSens[cellI] > 1.0) flameSens[cellI] = 1.0;
    }
    else if (XiPV[cellI] < 0.0) 
    {
         XiPV[cellI] = 0.0;
         flameSens[cellI] = 0.0; 
    }
    else if (XiPV[cellI] > 1.0) 
    {
         XiPV[cellI] = 1.0;
         flameSens[cellI] = 0.0; 
    }
    
    //- Calculate the thickening factor
    if (dynamicTh) 
    {
       TF[cellI] = 1.0 + (TFmax[cellI] - 1.0) * flameSens[cellI];
    }
    else 
    {
       TF[cellI] = TFmax[cellI];
    }
}

flameSens.correctBoundaryConditions();
TF.correctBoundaryConditions();

#include "evaluateLaplacian.H"
uPrime = 2.0 * mag(pow(LESdelta, 3) * fvc::curl(laplaceSim));

Info << "Calculating efficiency function for ATF" << endl;
efficiency->correct();
