/*------------------------------------------------------------------------------------*\
|                           Calculate Laplacian                                        |
|--------------------------------------------------------------------------------------|
| Velocity on LES mesh is mapped to a coarser evalLaplace mesh, the Laplacian operator |
| is evaluated on the coarse mesh and then inerpolated back to LES mesh                |
\*------------------------------------------------------------------------------------*/

//- Mapping of the velocity to a coarser evalLaplace mesh
forAll(cellsInLaplacianCell,cellLaplaceI) {
    scalar Umean[3];
    scalar rhoAverage = 0.0;
    for (label i = 0; i < 3; i++) {
        Umean[i] = 0.0;
    }
    forAll(cellsInLaplacianCell[cellLaplaceI],cellDNSI) {
        rhoAverage = rhoAverage + rho[cellsInLaplacianCell[cellLaplaceI][cellDNSI]];
    for (label i = 0; i < 3; i++) {
        Umean[i] = Umean[i] + rho[cellsInLaplacianCell[cellLaplaceI][cellDNSI]]*U[cellsInLaplacianCell[cellLaplaceI][cellDNSI]][i];
    }
    }
    // A coarse cell that no LES cell centre mapped into has rhoAverage == 0.
    // Dividing there wrote NaN into Uint, and the NaN then propagated through
    // coarseLap -> laplaceSim -> uPrime -> vb -> the mixing time scale.
    if (rhoAverage > SMALL) {
        for (label i = 0; i < 3; i++) {
            Uint[cellLaplaceI][i] = Umean[i] / rhoAverage;
        }
    } else {
        for (label i = 0; i < 3; i++) {
            Uint[cellLaplaceI][i] = 0.0;
        }
    }
}
Uint.correctBoundaryConditions();

//- The Laplacian is calculated on the coarser mesh
coarseLap = fvc::laplacian(Uint);
coarseLap.correctBoundaryConditions();

//- The Laplacian is interpolated back to LES mesh
forAll(cellsInLaplacianCell,cellLaplaceI) {
    forAll(cellsInLaplacianCell[cellLaplaceI],cellDNSI) {
        laplaceSim[cellsInLaplacianCell[cellLaplaceI][cellDNSI]] = coarseLap[cellLaplaceI];
    }
}
laplaceSim.correctBoundaryConditions();
