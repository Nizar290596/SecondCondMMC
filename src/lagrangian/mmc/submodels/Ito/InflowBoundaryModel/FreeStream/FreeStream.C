/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "FreeStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::FreeStream
(
    const dictionary& dict,
    CloudType& owner
)
:
    InflowBoundaryModel<CloudType>(owner),
    
    patches_(),
    
    Cm_(owner.deltaM()),
    
    particleFluxAccumulators_(),

    faceRandThreshold_()
{

    Info << "Construct FreeStream inflow model for Pope particles." << endl;

    // Identify which patches to use
    DynamicList<label> patches;

    forAll(owner.mesh().boundaryMesh(), p)
    {
        const polyPatch& patch = owner.mesh().boundaryMesh()[p];
        
        if (patch.type() == polyPatch::typeName)
        {
            patches.append(p);
        }
    }

    patches_.transfer(patches);

    particleFluxAccumulators_.setSize(patches_.size());

    faceRandThreshold_.setSize(patches_.size());

    StochasticLib1& rndGen = owner.rndGen();

    forAll(patches_, p)
    {
        const polyPatch& patch = owner.mesh().boundaryMesh()[patches_[p]];

        particleFluxAccumulators_[p] =  Field<scalar>(patch.size(), 0.0);

        faceRandThreshold_[p] =  Field<scalar>(patch.size());

        forAll(patch, f)
        {
            faceRandThreshold_[p][f] = rndGen.Random();
        }
    }
}


template<class CloudType>
Foam::FreeStream<CloudType>::FreeStream(const FreeStream<CloudType>& ibc)
:
    InflowBoundaryModel<CloudType>(ibc)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::~FreeStream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FreeStream<CloudType>::inflow()
{

    CloudType& cloud(this->owner());

    const fvMesh& mesh(cloud.mesh());
    
    const fvMesh& superMesh(cloud.pManager().superMesh());
    
    const List<label>& superCellForCell(cloud.pManager().superCellForCell());

    const scalar deltaT = cloud.time().deltaTValue();

    StochasticLib1& rndGen = cloud.rndGen();

    label particlesInserted = 0;

    const volScalarField::Boundary& boundaryRho
    (
        cloud.rho().boundaryField()
    );

    const volVectorField::Boundary& boundaryU
    (
        cloud.U().boundaryField()
    );

    forAll(patches_, p)
    {

        label patchI = patches_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        //- Add mass to the accumulators: negative face area dot product with 
        //- the velocity to point flux into the domain.Take a reference to the 
        //- particleFluxAccumulator for this patch
        Field<scalar>& pFA = particleFluxAccumulators_[p];
      
        //- Update pFA
        Field<scalar> pFAnew = (boundaryU[patchI] & -patch.faceAreas() )*boundaryRho[patchI]*deltaT;

        pFA += max(pFAnew, 0.0);
        
        //- Volume accumulator; used for determining number of particles to insert. 
        //- Cannot be negative else attempts to insert particles at outflow boundaries
        Field<scalar> pVA = pFA / boundaryRho[patchI];
        
        //- Loop over all faces as the outer loop to avoid calculating geometrical 
        //- properties multiple times.
        forAll(patch, f)
        {

            labelList faceVertices = patch[f];

            label globalFaceIndex = f + patch.start();
      
            label cellI = mesh.faceOwner()[globalFaceIndex];
  
            scalar fA = mag(patch.faceAreas()[f]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
                (
                    mesh,
                    globalFaceIndex,
                    cellI
                );

            //- Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //- Determine number of particles to inject
            scalar& faceAccumulator = pFA[f];

            const scalar& faceAccumulatorNew = pFAnew[f];

            if (faceAccumulator > 0 && faceAccumulatorNew > 0)
            {
                scalar particleMass = boundaryRho[patchI][f] 
                    * superMesh.V()[superCellForCell[cellI]] / cloud.pManager().Npc();

                scalar wt = particleMass / Cm_;

                scalar rI = faceAccumulator / particleMass;

                label nI = label(rI);

                if ( rI - nI > faceRandThreshold_[p][f] ) nI++;

                for (label I = 0; I < nI; I++)
                {

                    // Choose a triangle to insert on, based on their cumulative area
                    scalar triSelection = rndGen.Random();

                    // Selected triangle
                    label selectedTriI = -1;

                    forAll(cTriAFracs, tri)
                    {
                        selectedTriI = tri;

                        if (cTriAFracs[tri] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle, using the
                    // method from:
                    // Generating Random Points in Triangles
                    // by Greg Turk
                    // from "Graphics Gems", Academic Press, 1990
                    // http://tog.acm.org/GraphicsGems/gems/TriPoints.c
            
                    const tetIndices& faceTetIs = faceTets[selectedTriI];
                
                    Random tRand(rndGen.IRandom(0,1000000000));
          
                    point p = faceTetIs.faceTri(mesh).randomPoint(tRand);

                    // use tetPointRef corresponding to triangle and then call pointToBarycentric(p) !!!!!!!!

                    barycentric bary = faceTetIs.tet(mesh).pointToBarycentric(p);

                    cloud.addNewParticle
                    (
                        bary,
                        cellI,
                        globalFaceIndex,
                        faceTetIs.tetPt(),
                        wt,
                        particleMass,
                        patchI,
                        f,
                        false
                    );
                
                    particlesInserted++;

                
                }

                if (nI > 0)
                {
                    //- Mass remaining for next timestep
                    faceAccumulator -= nI*particleMass;
   
                    //- Reset random flux threshold for the face if a particle has been inserted
                    faceRandThreshold_[p][f] = rndGen.Random();
                }
            }
        }
    }
    
    reduce(particlesInserted,sumOp<label>());
  
    Info << "\t" << particlesInserted << " particles inserted at inflow boundaries." << endl;
}


// ************************************************************************* //
