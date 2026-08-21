fvMesh filterMesh
(
    IOobject
    (
        "filterMesh",
        filterTime.timeName(),
        filterTime,
        IOobject::MUST_READ
    )
);



