/*:00000000
 * intLength.C
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 *
 *      BUG:
 */

#include "intLength.H"

namespace Foam
{

intLength::intLength
(
    const volVectorField &U_,
    const fvMesh& Mesh_,
    statMesh& statProps_
)
:
    U(U_),
    Mesh(Mesh_),
    statProps(statProps_),
    intL(vector::zero)
{}

void intLength::updateL(){
    const label corLayers(statProps.corLayers());
    const label cellsPatch(statProps.cellsPatch());
    const label cellsCorrelation(statProps.cellsCorrelation());
    const label skipCells(statProps.skipCells());

    List<vector> Ravg(cellsCorrelation, vector::zero);    // autocorrelation in space Pope p.68 R=Ucorr/(Unorm1)^(1/2)
    List<vector> Ucorr(cellsCorrelation, vector::zero);    // <u(x)*u(x+r)>
    vector Unorm(vector::zero);                            // <u(x)*u(x)>
    vector intR(vector::zero);                            // vector to store integral length for all 3 velocity components
    scalar alpha(0),beta(0),meanL(0),n(0);                // needed for averaging
    label centerCell,neighbourcell,offset;                // index center cell = x, neighbour = x+i

    for(int layer=0; layer<corLayers; layer++){                    // iterate through all correlation layers
        offset = layer*skipCells*cellsPatch;                    // increase offset for the new slice
        for(int startFace=0;startFace<cellsPatch;startFace++){    // iterate through all faces on slice,
            centerCell    = startFace + offset;                    // set center cell on line

            n    += 1;
            alpha     = (n-1)/n;
            beta     = 1/n;

            Unorm[0] =alpha*Unorm[0]+beta*(U[centerCell][0]*U[centerCell][0]);
            Unorm[1] =alpha*Unorm[1]+beta*(U[centerCell][1]*U[centerCell][1]);
            Unorm[2] =alpha*Unorm[2]+beta*(U[centerCell][2]*U[centerCell][2]);

            for(int cell=0;cell<cellsCorrelation;cell++){                // cells along one line
                neighbourcell = centerCell + cell*cellsPatch;            // take center cell as start point
                Ucorr[cell][0]  = alpha * Ucorr[cell][0]  + beta * (U[centerCell][0]*U[neighbourcell][0]);
                Ucorr[cell][1]  = alpha * Ucorr[cell][1]  + beta * (U[centerCell][1]*U[neighbourcell][1]);
                Ucorr[cell][2]  = alpha * Ucorr[cell][2]  + beta * (U[centerCell][2]*U[neighbourcell][2]);
            }
        }
    }

    forAll(Ravg,cell){                                                    // calc the actual <Ravg>_slice
        for(int i=0;i<3;i++){                                             // unroll
            Ravg[cell][i] = Ucorr[cell][i]/(Unorm[i]);
        }
    }

    // Info << "Ravg" << Ravg << endl;

    for(int cell=1; cell < Ravg.size(); cell++){                   // integrate over <Ravg>
        intR = intR +  0.5*(Ravg[cell-1]+Ravg[cell])*statProps.deltaZ();
    }


    meanL = (intR[0]+intR[1]+intR[2])/3.0;
    Info << "integral length scale: " << intR <<  " mean: " << meanL <<  endl;


    if(intR[0]<0 || intR[1]<0  || intR[2]<0  ){
        Info << "negative lengthscales error ! Increasing the simulation time might help! Exiting ... " << endl;
    }

    intL = intR;        // store intL
}

scalar intLength::getMinL(){
    scalar minIntL(0);

    minIntL = min(min(intL[0],intL[1]),intL[2]);

    return minIntL;
}

scalar intLength::getMaxL(){
    scalar maxIntL(0);

    maxIntL = max(max(intL[0],intL[1]),intL[2]);

    return maxIntL;
}

scalar intLength::getMeanL(){
    scalar meanIntL(0);

    meanIntL = (intL[0]+intL[1]+intL[2])/3;

    return meanIntL;
}

}  // End namespace Foam
