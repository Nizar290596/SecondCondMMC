/*---------------------------------------------------------------------------*\
|                           Set enthalpy & species                            |
|-----------------------------------------------------------------------------|
| YEqvE and hEqvE are unique functions of the progress variable               |
| Set YEqvE and hEqvE based on the flamelet tables                            |
| Z42 mechanism: Larsson et al., Energy Fuels 31 (2017) 1904 – 1926           |
\*---------------------------------------------------------------------------*/

label inertIndex = -1;

volScalarField Yt1(0.0*YEqvE[0]);
volScalarField Temp(0.0*YEqvE[0]);
label FGMspecI = -1;
label FGMtempI = 1;
label pvRow;
scalar pvint0, pvint1;

//- Interpolate species mass fractions from the flamelet table (here use all species from Z42)
forAll(YEqvE, specieI)
{
    if (YEqvE[specieI].name() != inertSpecie)
    {
        volScalarField& Yi = YEqvE[specieI];
        Yi.max(0.0);
        if (YEqvE[specieI].name() == "O2")
        {
            FGMspecI = 3;
        }
        else if (YEqvE[specieI].name() == "H2O")
        {
            FGMspecI = 5;
        }
        else if (YEqvE[specieI].name() == "CH4")
        {
            FGMspecI = 4;
        }
        else if (YEqvE[specieI].name() == "CO2")
        {
            FGMspecI = 2;
        }
        else if (YEqvE[specieI].name() == "H2")
        {
            FGMspecI = 7;
        }
        else if (YEqvE[specieI].name() == "CO")
        {
            FGMspecI = 6;
        }
        else if (YEqvE[specieI].name() == "CH")
        {
            FGMspecI = 8;
        }
        else if (YEqvE[specieI].name() == "CH2")
        {
            FGMspecI = 9;
        }
        else if (YEqvE[specieI].name() == "CH2O")
        {
            FGMspecI = 10;
        }
        else if (YEqvE[specieI].name() == "CH3")
        {
            FGMspecI = 11;
        }
        else if (YEqvE[specieI].name() == "CH3O")
        {
            FGMspecI = 12;
        }
        else if (YEqvE[specieI].name() == "H")
        {
            FGMspecI = 13;
        }
        else if (YEqvE[specieI].name() == "H2O2")
        {
            FGMspecI = 14;
        }
        else if (YEqvE[specieI].name() == "HCO")
        {
            FGMspecI = 15;
        }
        else if (YEqvE[specieI].name() == "HO2")
        {
            FGMspecI = 16;
        }
        else if (YEqvE[specieI].name() == "NO")
        {
            FGMspecI = 17;
        }
        else if (YEqvE[specieI].name() == "O")
        {
            FGMspecI = 18;
        }
        else if (YEqvE[specieI].name() == "OH")
        {
            FGMspecI = 19;
        }

        forAll(YEqvE[specieI],cellI)
        {
            pvRow = floor(Xi.Vars("cLES").field()[cellI]*pvDisc);
            pvRow = min(pvRow,pvDisc);
            pvRow = max(0,pvRow);

            pvint0 = (pvRow + 0.0)/pvDisc;
            pvint1 = (pvRow + 1.0)/pvDisc;
 
            if (pvRow == pvDisc)
            {    
                Yi[cellI] = FGMtable[pvRow][FGMspecI];
            }
            else if (pvRow > 0 && pvRow < pvDisc)
            {    
                Yi[cellI] = interpolate2Table(Xi.Vars("cLES").field()[cellI], pvint0, pvint1, FGMtable[pvRow][FGMspecI], FGMtable[pvRow + 1][FGMspecI]);
            }
            else if (pvRow == 0)
            {
                Yi[cellI] = FGMtable[0][FGMspecI];
            }
        }
 
        Yi.max(0.0);
        Yt1 += Yi;
        
    } // loop over reacting species
    else if (YEqvE[specieI].name() == inertSpecie)
    {
        inertIndex = specieI;
    }
}

YEqvE[inertIndex] = scalar(1) - Yt1;
YEqvE[inertIndex].max(0.0);

forAll(YEqvE,specieI) 
{ 
    YEqvE[specieI].correctBoundaryConditions(); 
}

//- Interpolate temperature from the flamelet table
forAll(Temp,cellI)
{
    pvRow = floor(Xi.Vars("cLES").field()[cellI]*pvDisc);
    pvRow = min(pvRow,pvDisc);
    pvRow = max(0,pvRow);

    pvint0 = (pvRow + 0.0)/pvDisc;
    pvint1 = (pvRow + 1.0)/pvDisc;

    if (pvRow == pvDisc)
    {    
        Temp[cellI] = FGMtable[pvRow][FGMtempI];
    }
    else if (pvRow > 0 && pvRow < pvDisc)
    {    
        Temp[cellI] = interpolate2Table(Xi.Vars("cLES").field()[cellI], pvint0, pvint1, FGMtable[pvRow][FGMtempI], FGMtable[pvRow + 1][FGMtempI]);
    }
    else if (pvRow == 0)
    {
        Temp[cellI] = FGMtable[0][FGMtempI];
    }
}

//- Calculate the enthalpy based on the temperature from flamelet table
volScalarField hEqvEFGM = thermo.he(p,Temp);
hEqvE = hEqvEFGM;
volScalarField& T = thermo.T();

hEqvE.correctBoundaryConditions();
T.correctBoundaryConditions();
p.correctBoundaryConditions();
rho.correctBoundaryConditions();
thermo.correct(); 

rho = thermo.rho();
