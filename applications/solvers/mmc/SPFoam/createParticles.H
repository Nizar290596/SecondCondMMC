/*---------------------------------------------------------------------------*\
|                             Create Particles                                |
|-----------------------------------------------------------------------------|
| Creates the pope particle clouds                                            |
\*---------------------------------------------------------------------------*/

dictionary pdfMethod = mesh.solutionDict().subDict("PDFMETHOD");
Switch pdfMethodOn = pdfMethod.lookupOrDefault<Switch>("pdfMethodOn", false);
Switch iniReleaseOn = pdfMethod.lookupOrDefault<Switch>("iniReleaseOn", true);
label nSet = pdfMethod.lookupOrDefault<label>("nSet",1);

PtrList<basicReactingPopeCloud> pSets;
pSets.clear();
basicReactingPopeCloud* ptr=0;
        
if (!pdfMethodOn)
{
    Info << "PDF method is switched off. " << nl << endl;

}
else
{
    Info << "PDF method is switched on. Creating " << nSet << " set(s) of Pope particles. " << nl << endl;
    
    for(int i=1; i<=nSet; i++)
    {
        std::stringstream pName;
        pName<<"popeParticlesSet_"<<i;
    
        ptr = new basicReactingPopeCloud
        (
			pName.str(),
            mesh,
            U,
            DEff,
            rho,
            gradRho,
            Xi,
            pdfMethodOn,
            iniReleaseOn
        );
        pSets.append(ptr);
    }
    // Is typically already set in the constructor if readFields is called
    // Just calling it again for safety
    //basicReactingPopeCloud::particleType::setStaticProperties(pSets[0]);
    basicReactingPopeParticle::setStaticProperties(pSets[0]);
}
