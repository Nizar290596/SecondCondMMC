/*
 * intLength.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */


#ifndef INTLENGTH_H_
#define INTLENGTH_H_


#include "volFields.H"
#include "midPointSet.H"
#include "statMesh.H"

namespace Foam
{
    class intLength {
        private:

            const volVectorField& U;        // reference to velocity field
            const fvMesh &Mesh;             // auxMesh not really needed for uniform meshes 
            statMesh& statProps;
            vector intL;

        public:
            intLength(const volVectorField &U_, const fvMesh &Mesh_, statMesh& statProps_);

            void updateL();

            scalar getMinL();

            scalar getMeanL();

            scalar getMaxL();

            vector getL(){return intL;};
    };

}

#endif /* INTLENGTH_H_ */
