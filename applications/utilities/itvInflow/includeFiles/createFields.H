Info << "Reading field U from file itvInflowU";

autoPtr<volVectorField> U(new volVectorField
(
     IOobject
     (
         "itvInflowU",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh()
 ));
Info << "done ";

/*autoPtr<vectorField> origUField(
    new vectorField(
        U().boundaryField()[toBeCloned].internalField()
    )
*/
Info << "toBeCloned " <<  toBeCloned << endl;
fvPatchField<vector> origUField(U().boundaryField()[toBeCloned]);

//U.clear();                                    // free some memory

Info<< "creating field inflowU: \n" << endl;


wordList bTypes(6);

bTypes[0]="cyclic";        // TODO: implement fixedGradient !
bTypes[1]="cyclic";        // TODO: implement fixedGradient !
bTypes[2]="cyclic";
bTypes[3]="cyclic";
bTypes[4]="cyclic";
bTypes[5]="cyclic";

volVectorField inflowU
(
    IOobject
    (
        "inflowU",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    auxMesh,
    dimensionedVector("inflowU",dimensionSet(0,1,-1,0,0,0,0),vector::zero),
    bTypes
);

Info<< "creating field inflowP: \n" << endl;
volScalarField inflowP
(
    IOobject
    (
        "inflowP",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    auxMesh,
    dimensionedScalar("inflowP",dimensionSet(0,2,-2,0,0,0,0),0),
    bTypes
);

Info<< "creating field phi: \n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        auxMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(inflowU) & auxMesh.Sf()
);

Info<< "creating field nu: \n" << endl;
dimensionedScalar nu
(
   "nu",
   dimensionSet(0,2,-1,0,0,0,0),
   1
);

Info<< "creating field divU: \n" << endl;
volScalarField divU // field to check divergence
(
    IOobject
    (
        "divU",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    auxMesh,
    dimensionedScalar("divU",dimensionSet(0,0,-1,0,0,0,0),0),
    bTypes
);
