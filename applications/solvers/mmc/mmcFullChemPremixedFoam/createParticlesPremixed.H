/*---------------------------------------------------------------------------*\
|                             Create Particles                                |
|-----------------------------------------------------------------------------|
| Creates the pope particle clouds                                            |
\*---------------------------------------------------------------------------*/

dictionary pdfMethod = mesh.solutionDict().subDict("PDFMETHOD");
Switch pdfMethodOn = pdfMethod.lookupOrDefault<Switch>("pdfMethodOn", false);
Switch iniReleaseOn = pdfMethod.lookupOrDefault<Switch>("iniReleaseOn", true);
label nSet = pdfMethod.lookupOrDefault<label>("nSet",1);

DynamicList<basicPremixedReactingPopeCloud*> pSets;
pSets.clear();
basicPremixedReactingPopeCloud* ptr=0;


//  Enquire if pdfScatPlotOn is on and size of scatter plot field
Switch pdfScatPlotOn = pdfMethod.lookupOrDefault<Switch>("pdfScatPlotOn", false);
label pdfScatPlotSize = pdfMethod.lookupOrDefault<label>("pdfScatPlotSize", 0);
        
DynamicList<basicPremixedReactingPopeCloud*> pScatSets;
pScatSets.clear();
basicPremixedReactingPopeCloud* ptrScat=0;
std::stringstream pScatName;

// Set normalization value for the progress variable  (Updated for flamelet)
IOdictionary cloudPropertiesPrem
(
    IOobject
    (
        "cloudProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary progVarNorm(cloudPropertiesPrem.subDict("progVarNorm"));   
scalar YCO2eq(readScalar(progVarNorm.lookup("YCO2eq")));

if (!pdfMethodOn)
{
    Info << "PDF method is switched off. " << nl << endl;

}
else
{
    Info << "PDF method is switched on. Creating " << nSet << " set(s) of Pope particles. " << nl << endl;
    
    for(int i=1;i<=nSet;i++)
    {
        std::stringstream pName;
        pName<<"popeParticlesSet_"<<i;
    
        ptr = new basicPremixedReactingPopeCloud
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
    
    basicPremixedReactingPopeParticle::setStaticProperties(*pSets[0]);

    if (pdfScatPlotOn)
    {
        Info << nl << "PDF scatter plot switched on. Creating " << nSet << " set(s) of Pope particles for scatter plot output on-the-fly. " << nl << endl;
        
        Info << "Estimated size of scatter plot field is " << pdfScatPlotSize << endl;
        
        for(int j=1;j<=nSet;j++)
        {
            pScatName<<"popeParticlesScatPlotSet_"<<j;
            
            ptrScat = new basicPremixedReactingPopeCloud
            (
                pScatName.str(),
                mesh,
                U,
                DEff,
                rho,
                gradRho,
                Xi,
                pdfMethodOn,
                iniReleaseOn
            );
            
            pScatSets.append(ptrScat);
        }
    }
}
