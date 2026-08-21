/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "massFluxBasedPatchInjectionBase.H"

#include "SubField.H"
#include "triPointRef.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFluxBasedPatchInjectionBase::massFluxBasedPatchInjectionBase
(
    const polyMesh& mesh,
    const word& patchName,
    const scalar& maxRadius
)
:
    patchName_(patchName),
    maxRadius_(maxRadius),
    patchId_(mesh.boundaryMesh().findPatchID(patchName_)),
    patchArea_(0.0),
    patchNormal_(),
    cellOwners_(),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0)
{
    if (patchId_ < 0)
    {
        FatalErrorInFunction
            << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << mesh.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    updateMesh(mesh);
}


Foam::massFluxBasedPatchInjectionBase::massFluxBasedPatchInjectionBase(const massFluxBasedPatchInjectionBase& pib)
:
    patchName_(pib.patchName_),
    maxRadius_(pib.maxRadius_),
    patchId_(pib.patchId_),
    patchArea_(pib.patchArea_),
    patchNormal_(pib.patchNormal_),
    cellOwners_(pib.cellOwners_),
    triFace_(pib.triFace_),
    triToFace_(pib.triToFace_),
    triCumulativeMagSf_(pib.triCumulativeMagSf_),
    sumTriMagSf_(pib.sumTriMagSf_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massFluxBasedPatchInjectionBase::~massFluxBasedPatchInjectionBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFluxBasedPatchInjectionBase::updateMesh(const polyMesh& mesh)
{
    // Set/cache the injector cells
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
//    points_ = patch.points();
    const pointField& points = patch.points();

    cellOwners_ = patch.faceCells();
    Info << cellOwners_.size() << endl;
    


    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFace(2*patch.size());
    DynamicList<face> tris(5);
//    DynamicList<vector> patchNormal;
//    DynamicList<label> cellOwners;
//    DynamicList<point> patchPoints;

    // Set zero value at the start of the tri area list
    triMagSf.append(0.0);

    patchArea_ = 0.0;

    forAll(patch, facei)
    {    
        const face& f = patch[facei];
//        Info << mag(patch.faceCentres()[facei]) << endl;
        if(sqrt(pow(patch.faceCentres()[facei][0],2)+pow(patch.faceCentres()[facei][1],2))<=maxRadius_)
        {

            tris.clear();
            f.triangles(points, tris);

            forAll(tris, i)
            {
                triToFace.append(facei);
                triFace.append(tris[i]);
                triMagSf.append(tris[i].mag(points));
            }
//            patchPoints.append(patch.points()[facei]);
//            cellOwners.append(patch.faceCells()[facei]);
//           
	    scalar faceArea = mag(patch.faceAreas()[facei]);

            patchArea_ += faceArea;
//            patchNormal.append(patch.faceAreas()[facei]/faceArea);
        }
    }

    forAll(sumTriMagSf_, i)
    {
        sumTriMagSf_[i] = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::broadcast(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); i++)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // Convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); i++)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }

    const scalarField magSf(mag(patch.faceAreas()));

//    patchArea_ = sum(magSf);
    patchNormal_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());

    Info << "Injection area: " << patchArea_ << endl;
}


void Foam::massFluxBasedPatchInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    Random& rnd,
    volScalarField& massFlux,
    scalar& maxMassFlux,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    label hold=0;
    
    do
    {
        scalar areaFraction = rnd.globalPosition(scalar(0), patchArea_);

        if (cellOwners_.size() > 0)
        {
            // Determine which processor to inject from
            label proci = 0;
            forAllReverse(sumTriMagSf_, i)
            {
                if (areaFraction >= sumTriMagSf_[i])
                {
                    proci = i;
                    break;
                }
            }

            if (Pstream::myProcNo() == proci)
            {
                // Find corresponding decomposed face triangle
                label trii = 0;
                scalar offset = sumTriMagSf_[proci];
                forAllReverse(triCumulativeMagSf_, i)
                {
                    if (areaFraction > triCumulativeMagSf_[i] + offset)
                    {
                        trii = i;
                        break;
                    }
                }
                // Set cellOwner
                facei = triToFace_[trii];
                cellOwner = cellOwners_[facei];

                // Find random point in triangle
                const polyPatch& patch = mesh.boundaryMesh()[patchId_];
                const pointField& points = patch.points();
                const face& tf = triFace_[trii];
                const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
                const point pf(tri.randomPoint(rnd));

                // Position perturbed away from face (into domain)
                const scalar a = rnd.position(scalar(0.1), scalar(0.5));
                const vector& pc = mesh.cellCentres()[cellOwner];
                const vector d =
                    mag((pf - pc) & patchNormal_[facei])*patchNormal_[facei];
                scalar massFluxFaceI = massFlux.boundaryField()[patchId_][facei];
                position = pf - a*d;

                // Try to find tetFacei and tetPti in the current position
                mesh.findTetFacePt(cellOwner, position, tetFacei, tetPti);

                // tetFacei and tetPti not found, check if the cell has changed
                if (tetFacei == -1 ||tetPti == -1)
                {
                    mesh.findCellFacePt(position, cellOwner, tetFacei, tetPti);
                }
  
                // Both searches failed, choose a random position within
                // the original cell
                if (tetFacei == -1 ||tetPti == -1)
                {
                    // Reset cellOwner
                    cellOwner = cellOwners_[facei];
                    const scalarField& V = mesh.V();

                    // Construct cell tet indices
                    const List<tetIndices> cellTetIs =
                        polyMeshTetDecomposition::cellTetIndices(mesh, cellOwner);

                    // Construct cell tet volume fractions
                    scalarList cTetVFrac(cellTetIs.size(), 0.0);
                    for (label teti=1; teti<cellTetIs.size()-1; teti++)
                    {
                        cTetVFrac[teti] =
                            cTetVFrac[teti-1]
                          + cellTetIs[teti].tet(mesh).mag()/V[cellOwner];
                    }
                    cTetVFrac.last() = 1;

                    // Set new particle position
                    const scalar volFrac = rnd.sample01<scalar>();
                    label teti = 0;
                    forAll(cTetVFrac, vfI)
                    {
                       if (cTetVFrac[vfI] > volFrac)
                       {
                            teti = vfI;
                            break;
                        }
                    }

                    position = cellTetIs[teti].tet(mesh).randomPoint(rnd);
                    tetFacei = cellTetIs[teti].face();
                    tetPti = cellTetIs[teti].tetPt();
                }
                scalar rndNum = rnd.sample01<scalar>();

                if(rndNum*maxMassFlux<massFluxFaceI) //maxMassFlux_
                {
                    hold =1;
                }


            }
            else
            {
                cellOwner = -1;
                tetFacei = -1;
                tetPti = -1;

                // Dummy position
                position = pTraits<vector>::max;
            }
        }
        else
        {
            cellOwner = -1;
            tetFacei = -1;
            tetPti = -1;

            // Dummy position
            position = pTraits<vector>::max;
        }
        reduce(hold, sumOp<label>());
    } while (hold==0);
}


// ************************************************************************* //
