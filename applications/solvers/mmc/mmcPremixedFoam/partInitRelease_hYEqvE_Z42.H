/*---------------------------------------------------------------------------*\
|      Set enthalpy & species for initialization of MMC particles             |
|-----------------------------------------------------------------------------|
| Provides thin Eulerian fields for the initial release of the particles      |
| (the Eulerian fields for LES are made thick again in flamelet_hYEqvE_Z42.H) |
| Set YEqvE and hEqvE based on the flamelet tables                            |
| Z42 mechanism: Larsson et al., Energy Fuels 31 (2017) 1904 – 1926           |
\*---------------------------------------------------------------------------*/

Info << "The MMC particles' composition fields are initialized assuming thin particle flame structure." << nl << endl;

label inertIndexInit = -1;

volScalarField Yt1Init(0.0*YEqvE[0]);
volScalarField TempInit(0.0*YEqvE[0]);
label FGMspecIInit = -1;
label FGMtempIInit = 1;
label pvRowInit;
scalar pvint0Init, pvint1Init;

//- Interpolate species mass fractions from the flamelet table (here use all species from Z42)
forAll(YEqvE, specieI)
{
    if (YEqvE[specieI].name() != inertSpecie)
    {
        volScalarField& Yi = YEqvE[specieI];
        Yi.max(0.0);
        if (YEqvE[specieI].name() == "O2")
        {
            FGMspecIInit = 3;
        }
        else if (YEqvE[specieI].name() == "H2O")
        {
            FGMspecIInit = 5;
        }
        else if (YEqvE[specieI].name() == "CH4")
        {
            FGMspecIInit = 4;
        }
        else if (YEqvE[specieI].name() == "CO2")
        {
            FGMspecIInit = 2;
        }
        else if (YEqvE[specieI].name() == "H2")
        {
            FGMspecIInit = 7;
        }
        else if (YEqvE[specieI].name() == "CO")
        {
            FGMspecIInit = 6;
        }
        else if (YEqvE[specieI].name() == "CH")
        {
            FGMspecIInit = 8;
        }
        else if (YEqvE[specieI].name() == "CH2")
        {
            FGMspecIInit = 9;
        }
        else if (YEqvE[specieI].name() == "CH2O")
        {
            FGMspecIInit = 10;
        }
        else if (YEqvE[specieI].name() == "CH3")
        {
            FGMspecIInit = 11;
        }
        else if (YEqvE[specieI].name() == "CH3O")
        {
            FGMspecIInit = 12;
        }
        else if (YEqvE[specieI].name() == "H")
        {
            FGMspecIInit = 13;
        }
        else if (YEqvE[specieI].name() == "H2O2")
        {
            FGMspecIInit = 14;
        }
        else if (YEqvE[specieI].name() == "HCO")
        {
            FGMspecIInit = 15;
        }
        else if (YEqvE[specieI].name() == "HO2")
        {
            FGMspecIInit = 16;
        }
        else if (YEqvE[specieI].name() == "NO")
        {
            FGMspecIInit = 17;
        }
        else if (YEqvE[specieI].name() == "O")
        {
            FGMspecIInit = 18;
        }
        else if (YEqvE[specieI].name() == "OH")
        {
            FGMspecIInit = 19;
        }

        forAll(YEqvE[specieI],cellI)
        {
            scalar cLESATF = 0.0;
            scalar cLESLoc = Xi.Vars("cLES").field()[cellI];
            scalar TFmaxLoc = TFmax[cellI]; 
                    
            //- Rescaling of the reference space to particle flame thickness 
            if      (cLESLoc <= ((1 - 1/TFmaxLoc)/2))  { cLESATF = 0.0; }
            else if (cLESLoc >= ((1 + 1/TFmaxLoc)/2))  { cLESATF = 1.0; }
            else
            {
                scalar TFLoc = 0.0;
                         
                if (dynamicTh) { TFLoc = 1.0 + (TFmaxLoc - 1.0) * 16.0 * pow((cLESLoc * (1.0 - cLESLoc)), 2); }
                else           { TFLoc = TFmaxLoc; }
                         
                cLESATF = 0.5 + TFLoc * (cLESLoc - 0.5);
            }
            
            pvRowInit = floor(cLESATF*pvDisc);
            pvRowInit = min(pvRowInit,pvDisc);
            pvRowInit = max(0,pvRowInit);

            pvint0Init = (pvRowInit + 0.0)/pvDisc;
            pvint1Init = (pvRowInit + 1.0)/pvDisc;
 
            if (pvRowInit == pvDisc)
            {    
                Yi[cellI] = FGMtable[pvRowInit][FGMspecIInit];
            }
            else if (pvRowInit > 0 && pvRowInit < pvDisc)
            {    
                Yi[cellI] = interpolate2Table(cLESATF, pvint0Init, pvint1Init, FGMtable[pvRowInit][FGMspecIInit], FGMtable[pvRowInit + 1][FGMspecIInit]);
            }
            else if (pvRowInit == 0)
            {
                Yi[cellI] = FGMtable[0][FGMspecIInit];
            } 
        }
 
        Yi.max(0.0);
        Yt1Init += Yi;
        
    } // loop over reacting species
    else if (YEqvE[specieI].name() == inertSpecie)
    {
        inertIndexInit = specieI;
    }
}


YEqvE[inertIndexInit] = scalar(1) - Yt1Init;
YEqvE[inertIndexInit].max(0.0);

forAll(YEqvE,specieI) 
{ 
    YEqvE[specieI].correctBoundaryConditions(); 
    YEqvE[specieI].write();  
}


//- Interpolate temperature from the flamelet table
forAll(TempInit,cellI)
{
    scalar cLESATF = 0.0;
    scalar cLESLoc = Xi.Vars("cLES").field()[cellI];
    scalar TFmaxLoc = TFmax[cellI]; 
                    
    //- Rescaling of the reference space to particle flame thickness 
    if      (cLESLoc <= ((1 - 1/TFmaxLoc)/2))  { cLESATF = 0.0; }
    else if (cLESLoc >= ((1 + 1/TFmaxLoc)/2))  { cLESATF = 1.0; }
    else
    {
        scalar TFLoc = 0.0;
                         
        if (dynamicTh) { TFLoc = 1.0 + (TFmaxLoc - 1.0) * 16.0 * pow((cLESLoc * (1.0 - cLESLoc)), 2); }
        else           { TFLoc = TFmaxLoc; }
                         
        cLESATF = 0.5 + TFLoc * (cLESLoc - 0.5);
    }
    
    pvRowInit = floor(cLESATF*pvDisc);
    pvRowInit = min(pvRowInit,pvDisc);
    pvRowInit = max(0,pvRowInit);

    pvint0Init = (pvRowInit + 0.0)/pvDisc;
    pvint1Init = (pvRowInit + 1.0)/pvDisc;

    if (pvRowInit == pvDisc)
    {    
        TempInit[cellI] = FGMtable[pvRowInit][FGMtempIInit];
    }
    else if (pvRowInit > 0 && pvRowInit < pvDisc)
    {    
        TempInit[cellI] = interpolate2Table(cLESATF, pvint0Init, pvint1Init, FGMtable[pvRowInit][FGMtempIInit], FGMtable[pvRowInit + 1][FGMtempIInit]);
    }
    else if (pvRowInit == 0)
    {
        TempInit[cellI] = FGMtable[0][FGMtempIInit];
    }
}

//- Calculate the enthalpy based on the temperature from flamelet table
volScalarField hEqvEFGMInit = thermo.he(p,TempInit);
hEqvE = hEqvEFGMInit;
volScalarField& TInit = thermo.T();

hEqvE.correctBoundaryConditions();
TInit.correctBoundaryConditions();
p.correctBoundaryConditions();
rho.correctBoundaryConditions();
thermo.correct(); 

rho = thermo.rho();
