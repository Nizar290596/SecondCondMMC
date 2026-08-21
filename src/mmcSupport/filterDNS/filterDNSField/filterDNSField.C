/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "filterDNSField.H"
#include "tetOverlapVolume.H"
#include "polyMeshTetDecomposition.H"
#include "tetrahedron.H"
#include "tetPoints.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * 

void Foam::filterDNSField::buildCellToFilterLists()
{
    // First find all local filter cells and mark all cells that could not 
    // be attributed to a filter cell

    // Get cell center position of the DNS mesh
    const volVectorField& pos = mesh().C();

    cellToFilterCell_.resize(pos.size());
    cellToFilterCellVolRatio_.resize(pos.size());
    cellToFilterCellProcessorID_.resize(pos.size());

    nCellsInFilterCell_.resize(filterMesh_.C().size(),0);

    filterCellToCell_.resize(filterMesh_.C().size());
    filterCellToRemoteCells_.resize(filterMesh_.C().size());

    tetOverlapVolume overlapEngine;

    // List of local DNS cells marked for remote processing
    DynamicList<Tuple2<label,scalar>> DNSCellrequiresRemoteHandling;


    forAll(pos,celli)
    {
        // Get all filter cells that overlap with the DNS cell 
        labelList filterCellsInDNSCell = 
            overlapEngine.overlappingCells(filterMesh_,mesh_,celli);

        // Add all local found cells and add up their overlapping volumes
        scalar totalVolOverlap = 0;
        scalar totalVolOverlapRatio = 0;
        const scalar cellVolume = mesh_.V()[celli];
        for (label filterCellI : filterCellsInDNSCell)
        {
            // Calculate volume overlap
            scalar volOverlap = interVol(celli,filterCellI);
            scalar volOverlapRatio = volOverlap/cellVolume;

            // Only add cells with a volume overlap greater than 1%
            if (volOverlapRatio > 0.01)
            {
                cellToFilterCellVolRatio_[celli].append(volOverlapRatio);
                cellToFilterCell_[celli].append(filterCellI);
                cellToFilterCellProcessorID_[celli].append(Pstream::myProcNo());
                nCellsInFilterCell_[filterCellI] += 1;

                totalVolOverlap += volOverlap;
                totalVolOverlapRatio += volOverlapRatio;

                filterCellToCell_[filterCellI].append(celli);
            }
        }
        
        // If the volume overlap of the DNS cell to filterCell is not 100%
        // There must be remote filter cells
        if (totalVolOverlapRatio < (0.999))
        {
            Tuple2<label,scalar> tmp;
            tmp.first() = celli;
            tmp.second() = cellVolume - totalVolOverlap;
            DNSCellrequiresRemoteHandling.append(tmp);
        }
    }

    if (Pstream::parRun())
    {
        buildCellToFilterListsParallel(DNSCellrequiresRemoteHandling);
    }
}


void Foam::filterDNSField::buildCellToFilterListsParallel
(
    const DynamicList<Tuple2<label,scalar>>& DNSCellrequiringRemoteHandling
)
{
    Info << "\tStart parallel handling..." << endl;
    // Get cell center position of the DNS mesh
    const volVectorField& pos = mesh().C();

    // Get the bounding boxes of the filter mesh
    List<boundBox> processorFilterMeshBoundingBoxes = 
        getProcessorBoundingBoxes(filterMesh_);


    // We will stream to the buffer the boundBox of the cells that require
    // remote processing. This only works if the cells have hex character!
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Data format to send:
    // Send a Tuple2<Tuple2<label,scalar>,DynamicList<tetPoints>> containing 
    // First:  Cell ID and left volume to find in filtering cells
    // Second: A list of tetPoints of the cell
    // =========
    // Unfortunately we cannot directly stream into the buffer and read with 
    // eof() function. This causes an error with the readRaw() that the data 
    // is not read correctly. Hence, we package it first into lists... 
    // stupid but works

    typedef Tuple2<Tuple2<label,scalar>,DynamicList<tetPoints>> cellInfoType;

    List<DynamicList<cellInfoType>> cellInfoToSend(Pstream::nProcs());

    // Print a warning if the number of remote cells is large
    List<scalar> ratioOfRemoteCellHandlingPerProcessor(Pstream::nProcs());
    ratioOfRemoteCellHandlingPerProcessor[Pstream::myProcNo()] = 
        DNSCellrequiringRemoteHandling.size()/scalar(pos.size());
    Pstream::gatherList(ratioOfRemoteCellHandlingPerProcessor);

    if (Pstream::master())
    {
        forAll(ratioOfRemoteCellHandlingPerProcessor,procI)
        {
            if (ratioOfRemoteCellHandlingPerProcessor[procI] > 0.1)
                Info << "\tNumber of DNS cells requiring remote processing "
                     << "is larger than 10% of the total mesh size on processor "
                     << procI << ".\n"
                     << "\tThis makes parallel handling very slow.\n"
                     << "\tConsider decomposing your filterMesh in "
                     << "alignment of the DNS mesh" << nl << endl;        
        }
    }

    DynamicList<tetPoints> cellTetPoints;

    // Loop over all cells that require remote handling
    for (const Tuple2<label,scalar>& tmp : DNSCellrequiringRemoteHandling)
    {
        const label celli = tmp.first();
        const scalar volToFind = tmp.second();

        // Check if there is an overlap to any of the processors
        forAll(processorFilterMeshBoundingBoxes,procI)
        {
            if (procI == Pstream::myProcNo())
                continue;
            
            if (processorFilterMeshBoundingBoxes[procI].contains(pos[celli]))
            {
                cellInfoType cellInfo;
                cellInfo.first().first() = celli;
                cellInfo.first().second() = volToFind;

                // Decompose cell into tets
                List<tetIndices> cellTetIndices = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    celli
                );

                // Loop over all tets
                
                for (tetIndices& tetIndex : cellTetIndices)
                {
                    // First get a tetPointRef and then a tetPoints...
                    tetPointRef tetRef = tetIndex.tet(mesh_);
                    tetPoints tet
                    (
                        tetRef.a(),
                        tetRef.b(),
                        tetRef.c(),
                        tetRef.d()
                    );
                    cellTetPoints.append(std::move(tet));
                }
                cellInfo.second() = cellTetPoints;
                cellInfoToSend[procI].append(std::move(cellInfo));

                cellTetPoints.clear();
            }
        }
    }

    // Now send the lists to the other processors
    forAll(cellInfoToSend,procI)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UOPstream pStream(procI,pBufs);
        pStream << cellInfoToSend[procI];
    }

    // Free up the memory
    cellInfoToSend.clear();

    // Send buffers
    pBufs.finishedSends();

    // New buffer to send the information back
    PstreamBuffers pBufsToSendBack(Pstream::commsTypes::nonBlocking);


    // Check if the received DNS cellInfo can be found in the local filterMesh
    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        if (procI == Pstream::myProcNo())
            continue;
        // Data structure is:
        // First():  cellID label of DNS cell on local processor
        // Second(): Tuple2<> with first filter cell ID and second vol found
        DynamicList<Tuple2<label,Tuple2<label,scalar>>> cellInfoToSendBack;

        UIPstream pStream(procI, pBufs);
        DynamicList<cellInfoType> receivedCellInfo;
        pStream >> receivedCellInfo;

        // Loop over all cells 
        for (auto& cellInfo : receivedCellInfo)
        {
            const label DNSCellID = cellInfo.first().first();
            scalar volToFind = cellInfo.first().second();
            const DynamicList<tetPoints>& tets = cellInfo.second();
            
            // Create a treeBoundBox to find all potential filter cells
            treeBoundBox cellBB;
            for (const tetPoints& tet : tets)
            {
                cellBB.add(tet);
            }
            const indexedOctree<treeDataCell>& treeA = filterMesh_.cellTree();

            labelList overlappingCells = treeA.findBox(cellBB);

            // Loop over all tetPoints and find the nearest cell
            scalar totalVolAdded = 0;
            const scalar totalVolToFind = volToFind;
            for (label filterCellI : overlappingCells)
            {
                // Calculate the overlap with this cell
                scalar overlapVolume = overlap(filterCellI,filterMesh_,tets);

                // If the cell adds more than 1% add it
                if (overlapVolume/volToFind > 0.01)
                {
                    Tuple2<label,Tuple2<label,scalar>> cellFound;
                    cellFound.first() = DNSCellID;
                    cellFound.second().first() = filterCellI;
                    cellFound.second().second() = overlapVolume;
                    cellInfoToSendBack.append(std::move(cellFound));
                    volToFind = volToFind - overlapVolume;
                    totalVolAdded += overlapVolume;

                    // Mark that for this filter cell a value has been found
                    nCellsInFilterCell_[filterCellI] += 1;
                    Tuple2<label,Tuple2<label,scalar>> filterCellInfo;
                    filterCellInfo.first() = DNSCellID;
                    filterCellInfo.second().first() = procI;
                    filterCellInfo.second().second() = overlapVolume;
                    filterCellToRemoteCells_[filterCellI].append(std::move(filterCellInfo));
                }

                if (totalVolAdded/totalVolToFind > 0.99)
                    break;
            }
        }

        UOPstream pStreamToSendBack(procI,pBufsToSendBack);
        pStreamToSendBack << cellInfoToSendBack;
    }


    // Clear the old buffer
    pBufs.clear();

    // Send information back
    pBufsToSendBack.finishedSends();

    // Assemble information
    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UIPstream pStream(procI, pBufsToSendBack);
        DynamicList<Tuple2<label,Tuple2<label,scalar>>> receivedCellInfo;
        pStream >> receivedCellInfo;

        for (auto& cellInfo : receivedCellInfo)
        {
            const label celli = cellInfo.first();
            const label filterCelli = cellInfo.second().first();
            const scalar volOverlapRatio = cellInfo.second().second()/mesh_.V()[celli];

            cellToFilterCellVolRatio_[celli].append(volOverlapRatio);
            cellToFilterCell_[celli].append(filterCelli);
            cellToFilterCellProcessorID_[celli].append(procI);
        }
    }
}


Foam::scalar Foam::filterDNSField::interVol
(
    const label celli,
    const label filterCellI
) const
{
    tetOverlapVolume overlapEngine;
 
    // Note: avoid demand-driven construction of cellPoints
    // treeBoundBox bbFilterCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);
    const UList<label>& cellFaces = filterMesh_.cells()[filterCellI];
    treeBoundBox bbFilterCell(filterMesh_.points(), filterMesh_.faces()[cellFaces[0]]);
    for (label i = 1; i < cellFaces.size(); ++i)
    {
        bbFilterCell.add(filterMesh_.points(), filterMesh_.faces()[cellFaces[i]]);
    }
 
    scalar vol = overlapEngine.cellCellOverlapVolumeMinDecomp
    (
        mesh_,
        celli,
        filterMesh_,
        filterCellI,
        bbFilterCell
    );
 
    return vol;
}


Foam::List<Foam::boundBox> Foam::filterDNSField::getProcessorBoundingBoxes
(
    const fvMesh& mesh
) const
{
    // Get all processor bounding boxes for a quick check if a cell may be 
    // located on a processor 
    List<point> boundingBoxMinOfProcessors(Pstream::nProcs());
    List<point> boundingBoxMaxOfProcessors(Pstream::nProcs());

    boundingBoxMinOfProcessors[Pstream::myProcNo()] = mesh.bounds().min();
    boundingBoxMaxOfProcessors[Pstream::myProcNo()] = mesh.bounds().max();

    Pstream::gatherList(boundingBoxMinOfProcessors);
    Pstream::scatterList(boundingBoxMinOfProcessors);

    Pstream::gatherList(boundingBoxMaxOfProcessors);
    Pstream::scatterList(boundingBoxMaxOfProcessors);
    
    // Construct the processor bounding boxes
    List<boundBox> processorBoundingBoxes(Pstream::nProcs());
    forAll(boundingBoxMinOfProcessors,proci)
    {
        processorBoundingBoxes[proci] = boundBox
        (
            boundingBoxMinOfProcessors[proci],
            boundingBoxMaxOfProcessors[proci]
        );
    }
    return processorBoundingBoxes;
}


Foam::scalar Foam::filterDNSField::overlap
(
    const label celli,
    const fvMesh& mesh,
    const DynamicList<tetPoints>& tets
) const
{
    tetPointRef::sumVolOp iop;

    const cell& cFacesA = mesh.cells()[celli];
    const point& ccA = mesh.cellCentres()[celli];
  
    forAll(cFacesA, cFA)
    {
        label faceAI = cFacesA[cFA];
        const face& fA = mesh.faces()[faceAI];

        bool ownA = (mesh.faceOwner()[faceAI] == celli);
        label tetBasePtAI = 0;
        const point& tetBasePtA = mesh.points()[fA[tetBasePtAI]];
        for (label tetPtI = 1; tetPtI < fA.size() - 1; tetPtI++)
        {
            label facePtAI = (tetPtI + tetBasePtAI) % fA.size();
            label otherFacePtAI = fA.fcIndex(facePtAI);
            label pt0I = -1;
            label pt1I = -1;
            if (ownA)
            {
                pt0I = fA[facePtAI];
                pt1I = fA[otherFacePtAI];
            }
            else
            {
                pt0I = fA[otherFacePtAI];
                pt1I = fA[facePtAI];
            }
            
            const tetPoints tetA
            (
                ccA,
                tetBasePtA,
                mesh.points()[pt0I],
                mesh.points()[pt1I]
            );
            

            // Loop over tets of cellB
            for (const tetPoints& tetB : tets)
            {
                tetTetOverlap(tetA, tetB, iop);
            }
        }
    }
    

    if (iop.vol_ > 1E-30)
    {
       return iop.vol_;
    }

    return 0.0;   
}


void Foam::filterDNSField::tetTetOverlap
(
    const tetPoints& tetA,
    const tetPoints& tetB,
    tetPointRef::sumVolOp& insideOp
) const
{
    static tetPointRef::tetIntersectionList insideTets;
    label nInside = 0;
    static tetPointRef::tetIntersectionList cutInsideTets;
    label nCutInside = 0;
    
    tetPointRef::storeOp inside(insideTets, nInside);
    tetPointRef::storeOp cutInside(cutInsideTets, nCutInside);
    tetPointRef::dummyOp outside;
    
    // Precompute the tet face areas and exit early if any face area is
    // too small
    static FixedList<vector, 4> tetAFaceAreas;
    static FixedList<scalar, 4> tetAMag2FaceAreas;
    tetPointRef tetATet = tetA.tet();
    for (label facei = 0; facei < 4; ++facei)
    {
        tetAFaceAreas[facei] = -tetATet.tri(facei).areaNormal();
        tetAMag2FaceAreas[facei] = magSqr(tetAFaceAreas[facei]);
        if (tetAMag2FaceAreas[facei] < ROOTVSMALL)
        {
            return;
        }
    }
    
    static FixedList<vector, 4> tetBFaceAreas;
    static FixedList<scalar, 4> tetBMag2FaceAreas;
    tetPointRef tetBTet = tetB.tet();
    for (label facei = 0; facei < 4; ++facei)
    {
        tetBFaceAreas[facei] = -tetBTet.tri(facei).areaNormal();
        tetBMag2FaceAreas[facei] = magSqr(tetBFaceAreas[facei]);
        if (tetBMag2FaceAreas[facei] < ROOTVSMALL)
        {
            return;
        }
    }
    
    // Face 0
    {
        vector n = tetBFaceAreas[0]/Foam::sqrt(tetBMag2FaceAreas[0]);
        plane pl0(tetBTet.tri(0).a(), n, false);
        tetA.tet().sliceWithPlane(pl0, cutInside, outside);
        if (nCutInside == 0)
        {
            return;
        }
    }
    // Face 1
    {
        vector n = tetBFaceAreas[1]/Foam::sqrt(tetBMag2FaceAreas[1]);
        plane pl1(tetBTet.tri(1).a(), n, false);
        nInside = 0;
        for (label i = 0; i < nCutInside; i++)
        {
            const tetPointRef t = cutInsideTets[i].tet();
            t.sliceWithPlane(pl1, inside, outside);
        }
        if (nInside == 0)
        {
            return;
        }
    }
    // Face 2
    {
        vector n = tetBFaceAreas[2]/Foam::sqrt(tetBMag2FaceAreas[2]);
        plane pl2(tetBTet.tri(2).a(), n, false);
        nCutInside = 0;
        for (label i = 0; i < nInside; i++)
        {
            const tetPointRef t = insideTets[i].tet();
            t.sliceWithPlane(pl2, cutInside, outside);
        }
        if (nCutInside == 0)
        {
            return;
        }
    }
    // Face 3
    {
        vector n = tetBFaceAreas[3]/Foam::sqrt(tetBMag2FaceAreas[3]);
        plane pl3(tetBTet.tri(3).a(), n, false);
        for (label i = 0; i < nCutInside; i++)
        {
            const tetPointRef t = cutInsideTets[i].tet();
            t.sliceWithPlane(pl3, insideOp, outside);
        }
    }    
}


// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * * *

Foam::filterDNSField::filterDNSField
(
    const fvMesh& mesh,
    const fvMesh& filterMesh
)
:
    mesh_(mesh),
    filterMesh_(filterMesh)
{
    // Construct the filter addressing
    buildCellToFilterLists();
}

// ************************************************************************* //

