/*---------------------------------------------------------------------------*\
|                          Create evalLaplace Mesh                            |
|-----------------------------------------------------------------------------|
| Creates the evalLaplace mesh for the evaluation of the Laplacian            |
| Each evalLaplaceMesh side has to be two times larger than basic mesh side   |
\*---------------------------------------------------------------------------*/

Info << "Creating mesh for the evaluation of the Laplacian" << endl;

fvMesh evalLaplaceMesh
(
    IOobject
    (
        "evalLaplace",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

const label numCells = mesh.nCells();
const label numLaplacianCells = evalLaplaceMesh.nCells();

labelList laplacianCellForCell(numCells);
List<DynamicList<label> > cellsInLaplacianCell(numLaplacianCells);

if (numLaplacianCells != (numCells/8))
    FatalErrorIn ("createEvalLaplaceMesh.H")
    << "numLaplacianCells != numCells/8" << nl
    << abort(FatalError);
    
forAll(mesh.cells(),celli)
{
    point p = mesh.cellCentres()[celli];
        
    label laplacianTetFacei = 0;
    label laplacianTetPti = 0;
    label laplacianCelli = 0;
        
    evalLaplaceMesh.findCellFacePt(p,laplacianCelli,laplacianTetFacei,laplacianTetPti);
        
    laplacianCellForCell[celli] = laplacianCelli;

    const label celliiL = celli;
        
    cellsInLaplacianCell[laplacianCelli].append(celliiL);
}
