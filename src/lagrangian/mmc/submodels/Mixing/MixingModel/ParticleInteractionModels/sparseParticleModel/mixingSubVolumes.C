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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "mixingSubVolumes.H"
#include "OFstream.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include <queue>

Foam::mixingSubVolumes::mixingSubVolumes
(
    const fvMesh& mesh,
    const dictionary& dict,
    const scalar& ri,
    const scalar& Xii
)
:
    mesh_(mesh),
    bounds_(mesh.points(),false),
    nSubVolumeSets_
    (
        min(dict.lookupOrDefault<label>("nSubVolumeSets",1000),2.0*Pstream::nProcs())
    ),
    directProcessorNeighbors_(constructDirectNeighbors()),
    allProcessorBounds_(getProcessorBoundingBox()),
    ri_(ri),
    Xii_(Xii)
{
    if (Pstream::parRun())
        genSubVolumeSets();
}


void Foam::mixingSubVolumes::genSubVolumeSets()
{
    subVolumeSets_.reserve(nSubVolumeSets_);
    for (label k=0; k < nSubVolumeSets_; k++)
    {
        subVolumeSets_.append(constructSubVolumeSet());
    }
}


Foam::List<Foam::List<Foam::label>>
Foam::mixingSubVolumes::constructSubVolumeSet()
{
    DynamicList<List<label>> subVolumeSet;

    label startProcID = std::round
        (
            rnd_.globalSample01<scalar>()*(Pstream::nProcs()-1)
        );

    // Keep track of alread included processors
    HashSet<label> foundProcessors;

    DynamicList<label> newSubVolume;
    createSubVolume(newSubVolume,foundProcessors,startProcID);

    // Compact the sub-volume
    compactSubVolume(newSubVolume,foundProcessors);

    // This determines the first sub volume
    subVolumeSet.append(newSubVolume);

    // Find the next sub volume by exploring the neighborhood
    label k=0;
    while (foundProcessors.size() != Pstream::nProcs())
    {
        DynamicList<label> newSubVolume = exploreNeighborhood
        (
            foundProcessors,
            subVolumeSet[k++ % subVolumeSet.size()]
        );

        if (newSubVolume.size() > 0)
            subVolumeSet.append(newSubVolume);
    }

    return subVolumeSet;
}


void Foam::mixingSubVolumes::expandSubVolume
(
    DynamicList<label>& subVolume,
    HashSet<label>& foundProcessors
)
{
    // Create copy of the original subVolume
    const DynamicList<label> origSubVolume = subVolume;

    for (label procI : origSubVolume)
    {
        // Add direct neighbor if they are not part of another sub-volume
        for (label k : directProcessorNeighbors_[procI])
        {
            if (foundProcessors.insert(k))
                subVolume.append(k);
        }
    }
}


bool Foam::mixingSubVolumes::createSubVolume
(
    DynamicList<label>& newSubVolume,
    HashSet<label>& foundProcessors,
    const label procID
)
{
    if (foundProcessors.insert(procID))
        newSubVolume.append(procID);

    // Create first sub-volume to start with
    DynamicList<label> neighProcs(directProcessorNeighbors_[procID]);
    for (auto i : neighProcs)
    {
        if (foundProcessors.insert(i))
            newSubVolume.append(i);
    }

    label subVolumeSize = newSubVolume.size();

    if (subVolumeSize == 0)
        return false;

    while(!checkSubVolume(newSubVolume,1.1))
    {        
        // Loop over neighboring processors to expand sub volume
        expandSubVolume(newSubVolume,foundProcessors);

        // Break loop if no new processors were added to the sub volume
        if (newSubVolume.size() == subVolumeSize)
            break;
        else
            subVolumeSize = newSubVolume.size();
    }

    return checkSubVolume(newSubVolume);
} 



Foam::DynamicList<Foam::label> 
Foam::mixingSubVolumes::exploreNeighborhood
(
    HashSet<label>& foundProcessors,
    const List<label>& subVolumeInit
)
{
    // Loop over all processors of subVolumeInit and check if a valid
    // sub-volume can be created from their neighboring processors

    // procID with largest number of processors in the sub volume
    // --> required if no valid sub-volume is found
    label procIDWithLargestProcessorCount = -1;
    label maxProcsInSubVol = -1;

    for (label procI : subVolumeInit)
    {
        // Go over all neighbors as a starting point
    
        for (label neighProcI : directProcessorNeighbors_[procI])
    {
        HashSet<label> copyFoundProcs = foundProcessors;
        DynamicList<label> newSubVolume;
            if (createSubVolume(newSubVolume,copyFoundProcs,neighProcI))
        {
            foundProcessors = copyFoundProcs;

            compactSubVolume(newSubVolume,foundProcessors);
            return newSubVolume;
        }
        else
        {
            if (newSubVolume.size() > maxProcsInSubVol)
            {
                    procIDWithLargestProcessorCount = neighProcI;
                maxProcsInSubVol = newSubVolume.size();
            }
        }
    }
    }

    // If no valid sub-volume has been found, use the one with the largest
    // processor count
    DynamicList<label> newSubVolume;
    createSubVolume(newSubVolume,foundProcessors,procIDWithLargestProcessorCount);

    return newSubVolume;
}


void Foam::mixingSubVolumes::compactSubVolume
(
    DynamicList<label>& subVolume,
    HashSet<label>& foundProcessors
)
{
    // reduce first the found sub volume
    reduceSubVolume(subVolume,foundProcessors);

    // First make the sub-volume as cubic as possible
    makeCubic(subVolume,foundProcessors);
}


void Foam::mixingSubVolumes::makeCubic
(
    DynamicList<label>& subVolume,
    HashSet<label>& foundProcessors
)
{
    // Create a bounding box around the new sub-volume and try to fill the 
    // holes to make it a cube.


    // If sub-volume is of size 1, there is nothing to fill
    if (subVolume.size() == 1)
        return;

    // Generate a 3D grid with x,y,z coordinates and a 3D list of which cells
    // are already occupied

    // Store all min/max values of one dimension of a cube in one array 
    DynamicList<scalar> x;  // x coordinates of the new grid
    DynamicList<scalar> y;  // y coordinates of the new grid
    DynamicList<scalar> z;  // z coordinates of the new grid

    List3D<bool> cellInSubVolume;
    gen3DGrid(subVolume,x,y,z,cellInSubVolume);

    // Loop over 3D grid and check if a cell is empty - if yes try to find a 
    // processor box that fills this void
    for (label i=0; i < x.size()-1; i++)
    {
        for (label j=0; j < y.size()-1; j++)
        {
            for (label k=0; k < z.size()-1; k++)
            {
                if (!cellInSubVolume(i,j,k))
                {
                    // Create copy of original sub volume
                    DynamicList<label> origSubVolume = subVolume;

                    // Calculate cell center of that box
                    const point cellCenter
                    (
                        0.5*(x[i+1]+x[i]),
                        0.5*(y[j+1]+y[j]),
                        0.5*(z[k+1]+z[k])
                    );

                    // Loop over all not found processors and check if they 
                    // could fill the box
                    forAll(allProcessorBounds_,procI)
                    {
                        if (!foundProcessors.test(procI))
                        {
                            if 
                            (
                                allProcessorBounds_[procI].contains(cellCenter)
                             && directPath(subVolume,foundProcessors,procI)
                            )
                            {
                                subVolume.append(procI);
                                foundProcessors.insert(procI);
                            }
                        }
                    }
                }
            }
        }
    }
}


void Foam::mixingSubVolumes::reduceSubVolume
(
    DynamicList<label>& subVolume,
    HashSet<label>& foundProcessors
)
{
    // Calculate the span in each direction
    const vector span = subVolumeSpan(subVolume);
    
    // Find the direction with the maximum extend which can be reduced
    label maxDir = 0;
    scalar maxDist = 0;
    forAll(span,i)
    {
        if (span[i] > maxDist)
        {
            maxDist = span[i];
            maxDir = i;
        }
    }

    // Only attempt to reduce if the max dist is at least 20% larger than 
    // the required value
    if (maxDist*Xii_/ri_ > 1.2)
    {
        // Sort the processors of the sub-volume in the maxDir direction
        DynamicList<std::pair<scalar,label>> sortedProcessors;
        for (auto e : subVolume)
        {
            sortedProcessors.append
            (
                std::pair<scalar,label>
                (
                    allProcessorBounds_[e].centre()[maxDir],
                    e
                )
            );
        }
        std::sort(sortedProcessors.begin(),sortedProcessors.end());

        boundBox subVolBox;
        
        subVolume.clear();

        label k = 0;
        for (;k < sortedProcessors.size(); k++)
        {
            if (subVolBox.span()[maxDir]*Xii_/ri_ > 1.0)
                break;
            
            label procI = sortedProcessors[k].second;
            subVolBox.add(allProcessorBounds_[procI]);
            subVolume.append(procI);
        }

        for (;k < sortedProcessors.size(); k++)
            foundProcessors.unset(sortedProcessors[k].second);
    }
    
    return;
}


void Foam::mixingSubVolumes::gen3DGrid
(
    const DynamicList<label>& subVolume,
    DynamicList<scalar>& x,
    DynamicList<scalar>& y,
    DynamicList<scalar>& z,
    List3D<bool>& cellInVolume
)
{
    for (auto i : subVolume)
    {
        x.append(allProcessorBounds_[i].min().x());
        x.append(allProcessorBounds_[i].max().x());

        y.append(allProcessorBounds_[i].min().y());
        y.append(allProcessorBounds_[i].max().y());

        z.append(allProcessorBounds_[i].min().z());
        z.append(allProcessorBounds_[i].max().z());
    }
    // Sort the array
    std::sort(x.begin(),x.end());
    std::sort(y.begin(),y.end());
    std::sort(z.begin(),z.end());

    // Remove duplicate entries
    DynamicList<scalar> xSorted;  // x coordinates of the new grid
    DynamicList<scalar> ySorted;  // y coordinates of the new grid
    DynamicList<scalar> zSorted;  // z coordinates of the new grid

    scalar lastX = GREAT;
    scalar lastY = GREAT;
    scalar lastZ = GREAT;
    forAll(x,i)
    {
        if (x[i] != lastX)
            xSorted.append(x[i]);
        lastX = x[i];
        if (y[i] != lastY)
            ySorted.append(y[i]);
        lastY = y[i];
        if (z[i] != lastZ)
            zSorted.append(z[i]);
        lastZ = z[i];
    }

    x = xSorted;
    y = ySorted;
    z = zSorted;

    // 3D grid storing center point of cells and if they are inside or
    // outside the volume
    cellInVolume.resize
    (
        xSorted.size()-1, ySorted.size()-1, zSorted.size()-1, false
    );


    // Loop over all processors and mark all grid cells in the processor boundary
    for (label n : subVolume)
    {
        for (label i=0; i < xSorted.size()-1; i++)
        {
            for (label j=0; j < ySorted.size()-1; j++)
            {
                for (label k=0; k < zSorted.size()-1; k++)
                {
                    // Jump if this cell is already marked as in volume
                    if (cellInVolume(i,j,k))
                        continue;
                    
                    // Calculate cell center
                    const point cellCenter
                    (
                        0.5*(xSorted[i+1]+xSorted[i]),
                        0.5*(ySorted[j+1]+ySorted[j]),
                        0.5*(zSorted[k+1]+zSorted[k])
                    );

                    if (allProcessorBounds_[n].contains(cellCenter))
                        cellInVolume(i,j,k) = true;

                }
            }
        }
    }
}


Foam::scalar Foam::mixingSubVolumes::sphericity(const DynamicList<label>& subVolume)
{
    // Loop over all processors in the sub-volume and add each processor bound
    // box. Consider overlapping areas by splitting the bounding box into
    // new smaller cubes. 

    // Store all min/max values of one dimension of a cube in one array 
    DynamicList<scalar> x;  // x coordinates of the new grid
    DynamicList<scalar> y;  // y coordinates of the new grid
    DynamicList<scalar> z;  // z coordinates of the new grid

    List3D<bool> cellInVolume;
    gen3DGrid(subVolume,x,y,z,cellInVolume);

    // Find volume and area
    // Loop over each cell and check if state of cell true inside or false outside
    // changes over the face. If yes it is an outside face and is added to the 
    // area

    scalar V=0;
    scalar A=0;

    for (label i=0; i < x.size()-1; i++)
    {
        for (label j=0; j < y.size()-1; j++)
        {
            for (label k=0; k < z.size()-1; k++)
            {
                if (cellInVolume(i,j,k))
                {
                    // Calculate volume of cell 
                    V += 
                        (x[i+1]-x[i])
                      * (y[j+1]-y[j])
                      * (z[k+1]-z[k]);
                    

                    scalar xFacingArea = (y[j+1]-y[j])*(z[k+1]-z[k]);
                    scalar yFacingArea = (x[i+1]-x[i])*(z[k+1]-z[k]);
                    scalar zFacingArea = (x[i+1]-x[i])*(y[j+1]-y[j]);

                    // Check the faces if they are an outside face
                    if (i>0 && i < x.size()-2)
                    {
                        if (cellInVolume(i-1,j,k) != cellInVolume(i,j,k))
                            A += xFacingArea;
                        if (cellInVolume(i,j,k) != cellInVolume(i+1,j,k))
                            A += xFacingArea;
                    }
                    else
                    {
                        A += 2*xFacingArea;
                    }

                    if (j>0 && j < y.size()-2)
                    {
                        if (cellInVolume(i,j-1,k) != cellInVolume(i,j,k))
                            A += yFacingArea;
                        if (cellInVolume(i,j,k) != cellInVolume(i,j+1,k))
                            A += yFacingArea;
                    }
                    else
                    {
                        A += 2*yFacingArea;
                    }

                    if (k>0 && k < z.size()-2)
                    {
                        if (cellInVolume(i,j,k-1) != cellInVolume(i,j,k))
                            A += zFacingArea;
                        if (cellInVolume(i,j,k) != cellInVolume(i,j,k+1))
                            A += zFacingArea;
                    }
                    else
                    {
                        A += 2*zFacingArea;
                    }
                }
            }
        }
    }

    // Calculate sphericity
    return std::pow(M_PI,0.33)*std::pow(6*V,0.66)/A;
}


Foam::scalar Foam::mixingSubVolumes::compactness(const DynamicList<label>& subVolume)
{
    // [1] Based on work of E. Bribiesca, "A measure of compactness
    //     for 3D shapes", 1999
    //
    // Loop over all processors in the sub-volume and add each processor bound
    // box. Consider overlapping areas by splitting the bounding box into
    // new smaller cubes. 


    if (subVolume.size() == 1)
        return 1;

    // Store all min/max values of one dimension of a cube in one array 
    DynamicList<scalar> x;  // x coordinates of the new grid
    DynamicList<scalar> y;  // y coordinates of the new grid
    DynamicList<scalar> z;  // z coordinates of the new grid

    List3D<bool> cellInVolume;
    gen3DGrid(subVolume,x,y,z,cellInVolume);

    // Find all internal faces of the subvolume and calculate 
    // mean contact area

    // Internal faces (called contact face in [1])
    scalar A=0;

    // Number of internal faces
    scalar nInternalFaces = 0;

    // Number of cubeoids
    scalar n=0;

    for (label i=0; i < x.size()-1; i++)
    {
        for (label j=0; j < y.size()-1; j++)
        {
            for (label k=0; k < z.size()-1; k++)
            {
                if (cellInVolume(i,j,k))
                {
                    n++;
                    
                    scalar xFacingArea = (y[j+1]-y[j])*(z[k+1]-z[k]);
                    scalar yFacingArea = (x[i+1]-x[i])*(z[k+1]-z[k]);
                    scalar zFacingArea = (x[i+1]-x[i])*(y[j+1]-y[j]);

                    // Check the faces if they are an outside face
                    if (i>0 && i < x.size()-2)
                    {
                        if (cellInVolume(i-1,j,k) != cellInVolume(i,j,k))
                        {
                            A += xFacingArea;
                            nInternalFaces++;
                        }
                        if (cellInVolume(i,j,k) != cellInVolume(i+1,j,k))
                        {
                            A += xFacingArea;
                            nInternalFaces++;
                        }
                    }
                    else
                    {
                        A += 2*xFacingArea;
                        nInternalFaces +=2;
                    }

                    if (j>0 && j < y.size()-2)
                    {
                        if (cellInVolume(i,j-1,k) != cellInVolume(i,j,k))
                        {
                            A += yFacingArea;
                            nInternalFaces++;
                        }
                        if (cellInVolume(i,j,k) != cellInVolume(i,j+1,k))
                        {
                            A += yFacingArea;
                            nInternalFaces++;
                        }
                    }
                    else
                    {
                        A += 2*yFacingArea;
                        nInternalFaces += 2;
                    }

                    if (k>0 && k < z.size()-2)
                    {
                        if (cellInVolume(i,j,k-1) != cellInVolume(i,j,k))
                        {
                            A += zFacingArea;
                            nInternalFaces++;
                        }
                        if (cellInVolume(i,j,k) != cellInVolume(i,j,k+1))
                        {
                            A += zFacingArea;
                            nInternalFaces++;
                        }
                    }
                    else
                    {
                        A += 2*zFacingArea;
                        nInternalFaces += 2;
                    }
                }
            }
        }
    }

    // Calculate mean area of internal faces
    const scalar a = A/nInternalFaces;

    // Minimum possible contact area
    const scalar AcMin = a*(n-1);
    const scalar AcMax = 3*(a*n-a*std::pow(n,2/3));

    // Return compactness
    return (A-AcMin)/(AcMax-AcMin);
}


bool Foam::mixingSubVolumes::checkSubVolume
(
    const List<label>& subVolume,
    const scalar scalingFactor
)
{
    // Return true if one direction is larger than the required value.
    // Hence the k-d tree would split in this physical dimension

    // Calculate the extend of the sub-volume
    const vector span = subVolumeSpan(subVolume);

    for (scalar e : span)
    {
        if (e*Xii_/ri_ > scalingFactor)
            return true;        
    }

    return false;
}


Foam::vector 
Foam::mixingSubVolumes::subVolumeSpan(const List<label>& subVolume)
{
    if (subVolume.size() == 0)
        return vector();
    

    // Create empty bounding box
    boundBox subVolBox;

    for (auto i : subVolume)
    {
        subVolBox.add(allProcessorBounds_[i]);
    }

    return subVolBox.span();
}


Foam::List<Foam::boundBox> 
Foam::mixingSubVolumes::getProcessorBoundingBox()
{
    // get the span of the local processor mesh
    // this returns the extend in each physical dimension x,y,z
    List<boundBox> allProcessorBoundingBox(Pstream::nProcs());

    allProcessorBoundingBox[Pstream::myProcNo()] = bounds_;

    Pstream::gatherList(allProcessorBoundingBox);
    Pstream::scatterList(allProcessorBoundingBox);

    return allProcessorBoundingBox;
}


Foam::List<Foam::List<Foam::label>>
Foam::mixingSubVolumes::constructDirectNeighbors()
{
    // Find all direct neighbors of this processor:
    HashSet<label> directNeighborSet;
    const fvPatchList& patches = mesh_.boundary();
    for (auto& patch :  patches)
    {
        if 
        (
            isA<processorFvPatch>(patch)
         && !isA<processorCyclicFvPatch>(patch)
        )
            directNeighborSet.insert
            (
               refCast<const processorFvPatch>(patch).neighbProcNo()
            );
    }
    
    List<List<label>> allDirectNeighbors(Pstream::nProcs());

    DynamicList<label> directNeighbor(directNeighborSet.size());
    for (auto& e : directNeighborSet)
    {
        directNeighbor.append(e);
    }

    allDirectNeighbors[Pstream::myProcNo()] = directNeighbor;
    Pstream::gatherList(allDirectNeighbors);
    Pstream::scatterList(allDirectNeighbors);
    
    return allDirectNeighbors;
}


bool Foam::mixingSubVolumes::directPath
(
    const DynamicList<label>& subVolume,
    const HashSet<label>& foundProcessors,
    const label procI
)
{
    // Check if there is a direct path between the processor procI and one
    // of the processors in the sub-volume
    // Use Breadth first search
    
    HashSet<label> subVolumeSet(subVolume);
    HashSet<label> visitedProcessors;
    std::queue<label> processorsToExplore;
    processorsToExplore.push(procI);

    while (!processorsToExplore.empty())
    {
        label p = processorsToExplore.back();
        processorsToExplore.pop();
        visitedProcessors.insert(p);

        for (label pi : directProcessorNeighbors_[p])
        {
            // if direct path is possible return true
            if (subVolumeSet.test(pi))
                return true;

            if (!foundProcessors.test(pi) && !visitedProcessors.test(pi))
            {
                processorsToExplore.push(pi);
            }
        }
    }

    return false;
}



const Foam::List<Foam::List<Foam::List<Foam::label>>>&
Foam::mixingSubVolumes::subVolumeSets()
{
    return subVolumeSets_;
}


const Foam::List<Foam::List<Foam::label>>&
Foam::mixingSubVolumes::getSubVolumeSet()
{
    label subVolumeSetID = std::round
        (
            rnd_.globalSample01<scalar>()*(nSubVolumeSets_-1)
        );
    
    return subVolumeSets_[subVolumeSetID];
}


const Foam::List<Foam::label>&
Foam::mixingSubVolumes::getSubVolume(const label procID)
{
    auto& subVolumeSet = getSubVolumeSet();
    // Search for processor ID in subVolumeSet
    for (auto& subVolume : subVolumeSet)
    {
        for (auto& i : subVolume)
        {
            if (i == procID)
            {
                return subVolume;
            }
        }
    }

    FatalError 
        << "Could not find a sub-volume containing processor " << procID
        << exit(FatalError); 
    
    // Return dummy to avoid compiler warnings
    return subVolumeSets_[0][0];
}


void Foam::mixingSubVolumes::writeSubVolumeSets()
{
    if (Pstream::master())
    {
        label i=0;
        for (auto& subVolumeSet : subVolumeSets_)
        {
            OFstream of("subVolume-"+Foam::name(i++)+".dat");
            // as it is difficult to read OpenFOAM lists with python without
            // an explicit reader we wrtite the data here in a more user
            // friendly format
            for (auto& subVolume : subVolumeSet)
            {
                for (auto& i : subVolume)
                    of << i<<" ";
                of << endl;
            }
        }
    }
}


