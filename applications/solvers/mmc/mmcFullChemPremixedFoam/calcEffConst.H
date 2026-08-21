/*---------------------------------------------------------------------------*\
|                     Calculate the efficiency constant                       |
|-----------------------------------------------------------------------------|
| Calculates a constant based on the particle flame propagation speed,        |
| substitutes the reference variable relaxation                               |
\*---------------------------------------------------------------------------*/

if (!dynEffConst || !pdfMethodOn)
{
// Do nothing, keep the effConstant same as you read form the file
}
else
{
    forAll(pSets,pSi)
    { 			
        label stepCountOut = 0;
        scalar pWtSum = 0.0;
        scalar cPxpWt = 0.0;
        scalar avercP = 0.0;

        forAllIter(basicPremixedReactingPopeCloud, *pSets[pSi], iter)
        {
			if ((iter().progVarInt() >= 0.475) && (iter().progVarInt() <= 0.525))
			{
				cPxpWt += iter().wt() * iter().XiC()[0];
				pWtSum += iter().wt();
			}
			if (iter().rlxStep() > stepCountOut) stepCountOut = iter().rlxStep();
	    } 

        reduce(pWtSum, sumOp<scalar>());
        reduce(cPxpWt, sumOp<scalar>());

        avercP = cPxpWt / (pWtSum + VSMALL);

		if ((avercP > maxPartProgVar) && (stepCountOut >= minNumOfSteps))
        {
            effConst = effConst + effConstChangeRate;
            stepCountOut = 1;
        }
        if (avercP > crMaxPartProgVar)
        {
            effConst = effConst + 0.5;
            stepCountOut = 1;
        }
        else if ((avercP < minPartProgVar) && (stepCountOut > 10))
        {
            effConst = effConst - effConstChangeRate;
            stepCountOut = 1;
        }
        else  stepCountOut++;

        forAllIter(basicPremixedReactingPopeCloud, *pSets[pSi], iter)
        {
            iter().rlxStep() = stepCountOut;
        }

        Info << "Efficiency constant: " << effConst << token::TAB << " step: "<< stepCountOut << endl;
        
        using std::ofstream;
        ofstream dynEffConstOut;
        dynEffConstOut.open("dynEffConst.dat", ios_base::app);
        
        if(Pstream::myProcNo() == 0)  dynEffConstOut << "\r\n" << effConst << token::TAB << stepCountOut << token::TAB << avercP;
        
        dynEffConstOut.close();
    }
}
