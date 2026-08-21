/*
 * writeOut.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#ifndef WRITEOUTDATA_H_
#define WRITEOUTDATA_H_

#include "volFields.H"
#include "OFstream.H"
#include "inflowProperties.H" 

namespace Foam
{

class writeOut {

private:
    const volVectorField& U;        // ref to the generated field
    const fvMesh & mesh;            // ref to the aux mesh
    const Foam::Time& runTime;        // needed for object registry
    const inflowProperties& infProp;    // 
    const volVectorField origU;        // 

public:
    writeOut(const volVectorField& U_,const fvMesh & mesh_,const Foam::Time &runTime_,const inflowProperties& infProp_,const volVectorField origU);
    void exportData();            // writes out the interpolated and superimposed data

    vector interpolate(const vector position);// U at position x            


};

}
#endif /* RADIALDATA_H_ */
