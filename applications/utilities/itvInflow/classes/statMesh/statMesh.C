/*
 * statMesh.C
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#include "statMesh.H"

namespace Foam
{
    statMesh::statMesh(const fvMesh& auxMesh,inflowProperties& inflowProp)
    :
    clonePatchLabel_(0/*auxMesh.boundaryMesh().findPatchID(inflowProp.patchName())*/),
    cellsPatch_(auxMesh.boundaryMesh()[clonePatchLabel_].size()),
    nCells_(auxMesh.V().size()),
    nCellsZ_(nCells_/cellsPatch_),
    auxMeshLength_(auxMesh.bounds().maxDim()),        // TODO try span
    deltaZ_(auxMeshLength_/nCellsZ_),
    cellsCorrelation_(floor(5*inflowProp.intL()/deltaZ_)),
    skipCells_(cellsCorrelation_/2),
    //corLayers_(nCellsZ_/skipCells_) // correlation cell would exceed the cell index for the last couples of corlayers
    corLayers_(floor((nCellsZ_ - cellsCorrelation_)/skipCells_))
    {
        printInfo();
    }
    
    void statMesh::printInfo(){
        Info    <<  "Number of cells along the integration line: " << cellsCorrelation_
            << "\nNumber of independent corellation Layers: " <<  corLayers_
            << "\nNumber of cells on the copied patch: " <<  cellsPatch_
            <<  endl;

    }

}
