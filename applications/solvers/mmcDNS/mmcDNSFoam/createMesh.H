/*---------------------------------------------------------------------------*\
|                                Create Mesh                                  |
|-----------------------------------------------------------------------------|
| Creates the base mesh of the solver                                         |
| Note: superMesh is read by constructing the ItoPopeCloud                    |
\*---------------------------------------------------------------------------*/

Info << "Creating base mesh" << nl << endl;

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

