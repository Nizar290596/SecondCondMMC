if(pdfScatPlotOn)
{
    Info << "Sampling Particles" << endl;
 
    forAll(pSets,pSi)
    {
        const label cloudSize = pSets[pSi]->size();
            
        label scatCloudSize = pScatSets[pSi]->size();
        
        scalar sizeRatio = scalar(pdfScatPlotSize) / scalar(cloudSize);
            
        scalar dT = runTime.deltaT().value();
            
        scalar eT = runTime.endTime().value();
        
        scalar rT = dT / eT * sizeRatio;
        
        forAllIter(basicPremixedReactingPopeCloud, *pSets[pSi], iter)
        {
            scatCloudSize = pScatSets[pSi]->size();
            
            if (scatCloudSize >= pdfScatPlotSize)
            {
                WarningIn(" ") 
                    << endl 
                    << "Pope particle scatter plot cloud exceeds upper size bound" << endl
                    << "No more scatter data will be stored." << nl
                    << endl;
                
                break;
            }
            
            StochasticLib1& rndNum = pSets[pSi]->rndGen();
            
            if (rndNum.Random() < rT)
            {
                basicPremixedReactingPopeParticle& p = iter();
                
                pScatSets[pSi]->append(p.clone().ptr());
            }
        }
    }
}
