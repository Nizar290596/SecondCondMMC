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

Class
    Foam::mixingSubVolumes
 
Description
    Divide the total domain into sub-volumes by combining processor domains.
    Provides mixing volume for the particle mixing model. 

    General notation:

        sub-volume: a number of processors forming one mixing region

        set of sub-volumes: the total set of sub-volumes which combined fill
                            the domain space without any holes

SourceFiles
    mixingSubVolumes.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef mixingSubVolumes_H
#define mixingSubVolumes_H

#include "dictionary.H"
#include "fvMesh.H"
#include "boundBox.H"
#include "Pstream.H"
#include "Ostream.H"
#include "Random.H"
#include "List3D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declaration
class mixingSubVolumes;

/*---------------------------------------------------------------------------*\
                Class mixingSubVolumes Declaration
\*---------------------------------------------------------------------------*/

class mixingSubVolumes
{
    private:
            //- Const reference to the local mesh
            const fvMesh& mesh_;

            //- Store the bounding box of local processor
            //  The function mesh.bounds() cannot be used as it returns the 
            //  bounding box of the reconstructed mesh!
            const boundBox bounds_; 

            //- Number of sub-volume sets generated 
            const label nSubVolumeSets_;

            //- List of direct processor neighbors 
            const List<List<label>> directProcessorNeighbors_;

            //- Store each processor centre and span
            //  List of pairs for each processor with 
            //  their centre (first entry) and span (second entry)
            const List<boundBox> allProcessorBounds_;


            //- Normalisation factors
            const scalar ri_;
            const scalar Xii_;

            //- List of sub-volumes
            DynamicList<List<List<label>>> subVolumeSets_;

            //- Random number generator
            Random rnd_;

        // Private Member Functions

            //- Expand sub-volume by adding neighboring processors
            void expandSubVolume
            (
                DynamicList<label>& subVolume,
                HashSet<label>& foundProcessors
            );

            //- Create a sub-volume starting from a processor ID
            //  Returns true if created subVolume fullfills the checkSubVolume
            //  criteria
            bool createSubVolume
            (
                DynamicList<label>& newSubVolume,
                HashSet<label>& foundProcessors,
                const label procID
            );

            //- Explore the neighborhood of an existing subVolume to 
            //- create a new sub-volume adjacent to it 
            DynamicList<label> exploreNeighborhood
            (
                HashSet<label>& foundProcessors,
                const List<label>& subVolumeInit
            );

            //- Compact the sub-volume to create more spherical volumes
            void compactSubVolume
            (
                DynamicList<label>& subVolume,
                HashSet<label>& foundProcessors
            );

            //- Fill the sub-volume to create a cubic object
            void makeCubic
            (
                DynamicList<label>& subVolume,
                HashSet<label>& foundProcessors
            );

            //- Reduce the  sub-volume to smallest possible space
            void reduceSubVolume
            (
                DynamicList<label>& subVolume,
                HashSet<label>& foundProcessors
            );

            //- Check if a given sub-volume confirms with the set criterias
            //- for sub-volumes
            //  scalingFactor: Specify a scalar larger than 1 to extend the 
            //                 criteria
            bool checkSubVolume
            (
                const List<label>& subVolume,
                const scalar scalingFactor = 1.0
            );

            //- Check if there is a direct path between a processor and a 
            //  sub-volume. 
            //  Processors that are already found are not considered
            bool directPath
            (
                const DynamicList<label>& subVolume, 
                const HashSet<label>& foundProcessors,
                const label procI
            );

            //- Calculate the sphericity of a subVolume based on their 
            //- bounding boxes
            scalar sphericity(const DynamicList<label>& subVolume);

            //- Calculate the compactness of the sub volume
            //  Based on work of E. Bribiesca, "A measure of compactness
            //  for 3D shapes", 1999
            scalar compactness(const DynamicList<label>& subVolume);

            //- Get the physical size of all processor domains
            //  Returns a list of structs of type boundingBox
            List<boundBox> getProcessorBoundingBox();

            //- Construct list of directly neighboring processors
            //  These are processors with a direct boundary and which can be 
            //  found using the processor boundary condition
            List<List<label>> constructDirectNeighbors();

            //- Create a set of sub-volumes
            //  This is a random process and a different set of sub-volumes
            //  is found every time it is called
            List<List<label>> constructSubVolumeSet();

            //- Generate n sub-volumes sets 
            void genSubVolumeSets();

            //- Calculate the physical extend of a sub-volume
            //  Returns the maximum span in the x,y, and z direction
            vector subVolumeSpan(const List<label>& subVolume);


            //- Generate a 3D grid with cuboids and a list marking which 
            //- cubes are filled by the sub-volume
            //  Used to calculate sphericity or compactness 
            void gen3DGrid
            (
                const DynamicList<label>& subVolume,
                DynamicList<scalar>& x,
                DynamicList<scalar>& y,
                DynamicList<scalar>& z,
                List3D<bool>& cellInVolume
            );


    public:

        // Constructor

            //- Construct from dictionary
            mixingSubVolumes
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const scalar& ri,
                const scalar& Xii
            );

        // Member Functions

            //- Return list of sub-volumes generated
            //  Typically this function is not needed and use getSubVolume()
            //  to access one sub-volume configuration by random
            const List<List<List<label>>>& subVolumeSets();

            //- Return a random sub-volume set configuration
            const List<List<label>>& getSubVolumeSet();

            //- Return the sub-volume for the given processor ID
            const List<label>& getSubVolume(const label procID);

            //- Write out sub-volume sets for visualisation
            void writeSubVolumeSets();


};

} // End of namespace Foam

#endif

// ************************************************************************* //