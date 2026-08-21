/*---------------------------------------------------------------------------*\
|                             Create Particles                                |
|-----------------------------------------------------------------------------|
| Creates the pope particle clouds                                            |
\*---------------------------------------------------------------------------*/

PtrList<basicReactingPopeCloud> pSets;

if (!pdfMethodOn)
{
    Info << "PDF method is switched off. " << nl << endl;
}
else
{
    Info << "PDF method is switched on. Creating " << nSet << " set(s) of Pope particles. " << nl << endl;
    
    for(int i=1;i<=nSet;i++)
    {
        word pName = "popeParticlesSet_"+std::to_string(i);
    
        auto ptr = new basicReactingPopeCloud
        (
            pName,
            mmcCloudMeshRef.cref(),
            UFiltered,
            DEffFiltered,
            rhoFiltered,
            gradRhoFiltered,
            XiFiltered,
            iniReleaseOn,
            readCloudFields
        );

        pSets.append(ptr);

        // Is typically already set in the constructor if readFields is called
        // Just calling it again for safety
        basicReactingPopeCloud::particleType::setStaticProperties(pSets[0]);
    }
}
