/*
 * variances.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#ifndef VARIANCES_H_
#define VARIANCES_H_

#include "volFields.H"

namespace Foam
{

class variances {

private:
    volVectorField& U;                // reference to velocity field
    const fvMesh& Mesh;                // reference to mesh
    scalar MeshV;                    // the total mesh volume

    vector variance;

public:
    variances(volVectorField &U);

    void update();                            // updates the variance vector

    void normalise();                        // scales U field to obtain an unity variance

    void normalise2();                        // scales U field to obtain an unity variance

    void scale(tensor rms);                // applies lund's matrix transformation

    void corTensor();
};

}
#endif /* VARIANCES_H_ */
