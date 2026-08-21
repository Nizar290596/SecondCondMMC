/*
 * statMesh.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 *
 *    A simple class to hold inflow properties to
 *    assure data consistency
 */

#ifndef    STATMESH_H_
#define STATMESH_H_

#include "volFields.H"
#include "OFstream.H"
#include "inflowProperties.H"

namespace Foam
{

class statMesh {

private:

    const label clonePatchLabel_;// id of cloned patch
    const label cellsPatch_;    // needed ?
    const scalar nCells_;        // mesh size
    const scalar nCellsZ_;        // mesh size
    const scalar auxMeshLength_;
    const scalar deltaZ_;        //
    const label cellsCorrelation_;    // number of cells used for correlation along one line , needed ?
    const label skipCells_;        // number of cells between two corellation lines longitudinal
    const label corLayers_;        // number of corellation lines along the extruded mesh
    scalar lambda_;                // patch bounding box ratio

public:

    statMesh(const fvMesh& auxMesh,inflowProperties& inflowProp);

    void printInfo();

    label cellsPatch(){return cellsPatch_;}      

    label cellsCorrelation(){return cellsCorrelation_;}

    label skipCells(){return skipCells_;}           

    label corLayers(){return corLayers_;}

    scalar deltaZ(){return deltaZ_;}

    label nCellsZ(){return nCellsZ_;}

    void setLambda(scalar& lambda){lambda_ = lambda;}

    scalar getLambda(){return lambda_;}
};
}
#endif /* RADIALDATA_H_ */
