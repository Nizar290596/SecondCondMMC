/*---------------------------------------------------------------------------*\
|                             Flamelet Closure                                |
|-----------------------------------------------------------------------------|
| Updates the progress variable source term based on the flamelet tables      |
\*---------------------------------------------------------------------------*/

forAll(pvSource,cellI)
{
    label pvRow;
    scalar pvint0, pvint1;

    pvRow = floor(Xi.Vars("cLES").field()[cellI]*pvDisc);
    pvRow = min(pvRow,pvDisc);
    pvRow = max(0,pvRow);

    pvint0 = (pvRow + 0.0)/pvDisc;
    pvint1 = (pvRow + 1.0)/pvDisc;

    if (pvRow == pvDisc)
    {    
        pvSource[cellI] = FGMtable[pvRow][sourceI];
    }
    else if (pvRow > 0 && pvRow < pvDisc)
    {    
        pvSource[cellI] = interpolate2Table(Xi.Vars("cLES").field()[cellI], pvint0, pvint1, FGMtable[pvRow][sourceI], FGMtable[pvRow + 1][sourceI]);
    }
    else 
    {
        pvSource[cellI] = 0.0;
    }
}

pvSource.correctBoundaryConditions();
