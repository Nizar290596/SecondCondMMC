/*---------------------------------------------------------------------------*\
Todo:    1. set boundary conditions for inflowU
        2. save exMesh in systems folder instead in timestep folder
        3.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "boundaryMesh.H"
#include "meshSearch.H"
#include <time.h>
#include "classes/intLength/intLength.H"
#include "classes/variances/variances.H"
#include "classes/radialData.H"
#include "statMesh.H"
#include "inflowProperties.H"
#include "interpolation.H"
#include "includeFiles/writeFiles.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("nr", "don't randomize");
    argList::validOptions.insert("wI", "write out complete inflow data");
    

    // Part 1 Read Meshes/Data/Fields/Dictionaries
    #include "setRootCase.H"
    #include "createTime.H"
    #include "readInflowDict.H"                            // reads inflow dict !!
    inflowProperties inflowProp(runTime,inflowDict);    // create object holding inflow props
    #include "readMeshes.H"
    #include "createFields.H"
    statMesh statisticProps(auxMesh,inflowProp);        // create object holding info for statistics
    intLength intL(inflowU,auxMesh,statisticProps);        // initialise integral lentgh scale
    
    scalar timeOffset(0);
    bool restart(true);

    while(restart){
        restart = false;
        #include "randomizeField.H"                        // randomizes the inflowU field
        variances var(inflowU);                            // initialise variances
        var.normalise();

        // Info << "Diffusion Time Step:" << diffDeltaT << endl;
        runTime.setDeltaT(inflowProp.deltaT());

        // Part 2 Diffusion
        Info<< nl << "Starting Inflow Generation" << endl;

        volVectorField storedInflowU(inflowU);            // vector field to store converged components
        vector stored(vector::zero);
        intL.updateL();

        while(intL.getMinL() < inflowProp.intL())        // run until smallest scale has converged
        {
            #include "UEqn.H"

            if(args.found("wI"))
            {
                divU.write();
                runTime.write();
                inflowU.write();
            }
            runTime++;
            intL.updateL();

            for(int i=0;i<3;i++){
                if(stored[i]==0 &&  intL.getL()[i] > inflowProp.intL()){                                    // store converged
                    forAll(inflowU,cellI){
                        storedInflowU[cellI][i] = inflowU[cellI][i];
                    }
                }
            }
        }

        Info << "elapsed diffusion time: "
            <<  runTime.timeOutputValue() - runTime.startTime().value()
            << endl;

        var.corTensor();
        var.scale(inflowProp.var());                // apply Lund matrix transformation to match R-tensor
        var.corTensor();
        divU=fvc::div(inflowU);


        if(args.found("wI"))                    // write out last time step after transformation
        {
            divU.write();
            runTime.write();
            inflowU.write();
        }

         #include "writeOutData.H"
    }
}

// ************************************************************************* //
