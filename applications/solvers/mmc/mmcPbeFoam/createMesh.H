    Info << "Creating meshes" << nl << endl;

// **********************************************************************
//  -region defaultRegion
// **********************************************************************
    Info << token::TAB << "Creating mesh for fv solution" << nl << endl;
    fvMesh mesh
    (
        IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );


// **********************************************************************
//  -region super
// **********************************************************************
    Info << token::TAB << "Creating super mesh for Pope particle number control and density feedback" << nl << endl;
    fvMesh superMesh
    (
        IOobject
        (
            "super",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

// **********************************************************************
//  - mesh mapping
// **********************************************************************
    Info << token::TAB << "Creating mesh to superMesh correspondence " << nl << endl;
    

    Info << token::TAB << "Creating superCellForCell list " << nl <<endl;
    const label numCells = mesh.nCells();
    labelList superCellForCell(numCells);
    
    Info << token::TAB << "Creating cellsInSuperCell list " << endl;
    const label numSuperCells = superMesh.nCells();
    List<DynamicList<label> > cellsInSuperCell(numSuperCells);
    
    if (numSuperCells > numCells)
        FatalErrorIn
        (
            "createMesh.H"
        )
        << "numSuperCells > numCells" << nl
        << abort(FatalError);
    

    forAll(mesh.cells(),celli)
    {
        point p = mesh.cellCentres()[celli];
        
        label superTetFacei = 0;
        label superTetPti = 0;
        label superCelli = 0;
        
        superMesh.findCellFacePt(p,superCelli,superTetFacei,superTetPti);
        
        superCellForCell[celli] = superCelli;

        const label cellii = celli;
        
        cellsInSuperCell[superCelli].append(cellii);
    }
    
