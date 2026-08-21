// Read the settings of the PDF method from:
// system/fvSolution PDFMETHOD{} dictionary

dictionary pdfMethod = filterMesh.solutionDict().subDict("PDFMETHOD");
Switch pdfMethodOn = pdfMethod.lookupOrAddDefault<Switch>("pdfMethodOn", false);
Switch iniReleaseOn = pdfMethod.lookupOrAddDefault<Switch>("iniReleaseOn", true);
Switch readCloudFields = pdfMethod.lookupOrAddDefault<Switch>("readFields", true);
label nSet = pdfMethod.lookupOrAddDefault<label>("nSet",1);
Switch filterDNSFields = pdfMethod.lookupOrAddDefault<Switch>("filterDNSFields", true);

Info << "PDF Settings" << pdfMethod << endl;

// =============================================================================
// Note: If filterDNSField is switched off, the DNS velocity, density,
//       temperature etc. fields are used for the mmc cloud. This means, 
//       that the cloud must be constructed with the DNS mesh, to allow for 
//       the correct interpolation methods of the tracking data!
//       To allow the switch between DNS and LES mesh for the cloud, a new
//       reference pointer, the mmcCloudMesh is introduced. 
//       Depending on the value of the switch filterDNSField, the mmcCloudMesh
//       is either the DNS or filter/LES mesh.

refPtr<fvMesh> mmcCloudMeshRef;


if (filterDNSFields)
{
    mmcCloudMeshRef.cref(filterMesh);
}
else
{
    mmcCloudMeshRef.cref(mesh);
}




