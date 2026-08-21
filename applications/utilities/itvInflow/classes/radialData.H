/*
 * radialData.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#ifndef RADIALDATA_H_
#define RADIALDATA_H_

#include "volFields.H"
#include "OFstream.H"


namespace Foam
{

class radialData {

private:
    const volVectorField& U;                // reference to velocity field
    const fvMesh & mesh;                    // reference to mesh
    const Foam::Time& runTime;

public:
    radialData(const volVectorField& U_,const fvMesh & mesh_,const Foam::Time &runTime_,const polyPatch& patch,const label& cells);//

    fvPatchField<tensor> redToPatch(const polyPatch& patch,label cells);

    void redToLine(fvPatchField<tensor> varPatch);
    void writeOutVarPatch(const fvPatchField<tensor>& patch);
    void writeOutMeanPatch(const fvPatchField<vector>& patch);
};

}
#endif /* RADIALDATA_H_ */
