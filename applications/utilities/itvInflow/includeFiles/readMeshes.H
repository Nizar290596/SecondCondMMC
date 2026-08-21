Info << "Create mesh for time = " << runTime.timeName() << endl;

Foam::IOobject meshIO
(
        Foam::fvMesh::defaultRegion,
        runTime.constant(),
        runTime,
        Foam::IOobject::MUST_READ
);

autoPtr<fvMesh> mesh(new fvMesh(meshIO));


const polyBoundaryMesh& patches = mesh().boundaryMesh();        // a reference to the patches of the original mesh
const label toBeCloned = patches.findPatchID(inflowProp.patchName());    // find the patch which is to be cloned
const polyPatch& pp = patches[toBeCloned];                // set a reference to that patch



Foam::IOobject auxMeshIO
(
                "exMesh",
                runTime.constant(),
                runTime,
                Foam::IOobject::MUST_READ
);

Info << "reading extruded Mesh: ";
fvMesh auxMesh(auxMeshIO);
Info << " done " << endl;

    scalar nFaceJet = pp.size();
    Info << "face number of the jet patch is: " << nFaceJet << endl;

//find the patch of the xoy plane of itvDomain 
    const polyBoundaryMesh& itvpatches = auxMesh.boundaryMesh();

    const label xoy = itvpatches.findPatchID("patch0_half0"); 

    const polyPatch& itvpp = itvpatches[xoy];

    scalar nFaceItv = itvpp.size();
    Info << "face number of the itv patch is: " << nFaceItv << endl;

   // Info << "if statement: " << nFaceJet > (nFaceItv * 0.5) << endl;

//    if(nFaceJet > (nFaceItv * 0.5))
//  {      
//            FatalErrorIn
//            (
//                "  "
//            )   
//                << "Resolution of the exMesh if not sufficient." << nl
//                << "Generate exMesh in constant/ folder with dimesion (n n 6000), where n*n is"<<nl
//                << "suggested at least twice larger than the faces of jet patch. " 
//                << exit(FatalError);
//        }


