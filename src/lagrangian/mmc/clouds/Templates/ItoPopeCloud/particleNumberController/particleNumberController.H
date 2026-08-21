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
    Foam::particleNumberController
 
Description
    Controls the particle number density in mmc using a super mesh
    
    The super mesh is read from the global constant/ folder and does not need
    to be decomposed
    
    Also initializes a particle field by returning a list of positions and
    particle mass 
    

SourceFiles
    particleNumberController.C
    particleNumberControllerTemplates.H

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#ifndef particleNumberController_H
#define particleNumberController_H

#include "List.H"
#include "fvMesh.H"
#include "StochasticLib.H"
#include "fvCFD.H"
#include "barycentric.H"
#include "tetIndices.H"
#include "polyMeshTetDecomposition.H"
#include <algorithm>                    // for std::random_shuffle

namespace Foam
{

class particleNumberController 
{
    public:
    
        // Container for particle properties to initialize the particle field
        // Required for initializing particle population based on a field 
        struct pContainer
        {
            barycentric location;   // location in the cell
            label cell;             // local mesh cell ID
            label tetFace;          // local mesh tetFace
            label tetPt;            // local mesh tetPoint
            scalar mass;            // particle mass
        };

    private:

        // Private Member Variables

            //- Store reference to local mesh
            const fvMesh& mesh_;

            //- Particle management settings
            const dictionary particleManagementDict_;
            
            //- Switch on numbers control
            const Switch numCtrlOn_;
            
            //- Wanted numbers of particles in each super cell
            const label Npc_;
            
            //- Lower limit of particles in each super cell
            const label Nlo_;
            
            //- Upper limit of particles in each super cell
            const label Nhi_;
            
            
            //- Number of cells in the super mesh
            label nSuperCells_;
            
            //- Number of super cells on the processor
            label nSuperCellsOnProcessor_;
            
            //- Time instance pointing to the root case not the processor
            //  case
            autoPtr<Time> time_;
        
            //- super mesh
            autoPtr<fvMesh> superMesh_;
            
            //- Connection of cells to super cell
            List<label> superCellForCell_;
            
            //- Cells that are located in each super cell
            List<DynamicList<label> > cellsInSuperCell_;
            
            //- Is super cell located on current processor
            List<bool> superCellOnProcessor_;
            
            //- Processor ID for each super cell 
            //  First index is super cell, second index is processor
            List<List<label>> superCellToProcID_;            
        
        // Private Member Functions
    
            //- Construct mesh from global constant/ directory
            fvMesh* constructSuperMesh(const fvMesh& mesh);
            
            //- Generate the associations between super cell and processors
            //  A super cell can be split on multiple processors
            void constructSuperCellToProcID();
    
            //- Collect the cell positions contained in each super cell 
            //  Handles parallel cases with Pstream 
            List<DynamicList<point>> getCellPositionsInSuperCell();
            
            //- Synchronize random seed for all super cells
            //  This is required to ensure, that the same super cell uses the 
            //  same random number, even if it is distributed over multiple 
            //  processors
            List<int> synchronizeRandomSeed() const;
            
            //- Generate a random points within a super cell
            List<point> randomPointsInSuperCell
            (
                const label superCellI,
                const label nPoints,
                const int seed
            ) const;
            
            //- Convert a point into barycentric coordinates of a given cell
            barycentric transformPoint(const point& p, const tetPointRef tet);
    
        // Template functions
        
            //- scalar comparison function used for sorting 
            //- (should be mixture fraction but used mass for testing)
            template<class ParticleType>
            class lessWt
            {
                public:

                    inline bool operator()
                    (
                        Tuple2<scalar,ParticleType*>& one,
                        Tuple2<scalar,ParticleType*>& two
                    ) const {return one.first() < two.first();}
            };
            
            //- Assemble the particle list from neighboring processors if 
            //- necessary
            template<class CloudType>
            List
            <
                DynamicList
                <
                    Tuple2
                    <
                        scalar,
                        typename CloudType::particleType*    
                    >
                >
            > 
            collectParticleList(CloudType& cloud);
    
    public:
        // Constructors
        
            //- Default constructor
            //  Reads the superMesh from the constant folder
            particleNumberController
            (
                const fvMesh& mesh, 
                const dictionary& particleManagementDict
            );
            
            
        // Public Member functions 
        
            //- Correct particle number in each super cell
            //  Works only with ItoPopeParticle of the ItoPopeCloud as it 
            //  requires a weight field of the particle
            template<class CloudType>
            void correct
            (
                CloudType& cloud
            );

            //- Return list of particle pointers for each superCell
            //- -- Only considers particles in the cloud given
            template<class CloudType>
            List<DynamicList<typename CloudType::particleType*>>
            getParticlesInSuperCellList
            (
                CloudType& cloud
            ) const;

            //- Return list of particle pointers to the particles located in
            //- the superCellI -- Only considers particles in the cloud given
            template<class CloudType>
            DynamicList<typename CloudType::particleType*>
            getParticlesInSuperCell
            (
                CloudType& cloud, 
                const label superCellI
            ) const;

            //- Reconstructs for a given field a list of all values in each 
            //- super cell
            //  Access is List[superCellI][cellInSuperCell]
            //  Size of List[superCellI] can be larger than local mesh for 
            //  parallel cases
            template<class FieldType>
            List<DynamicList<typename FieldType::value_type>>
            fieldDataInSuperCell(const FieldType& field) const;

            //- Create a list of particle positions such that the total 
            //  mass of the particles matches the mass in the domain
            List<pContainer> initParticlePosBasedOnField
            (
                const volScalarField& rho
            );

            //- Return the super cell for the given cell index
            label getSuperCellID(const label cellID) const;

        // Inline functions 
        
            //- Access to superCellForCell list
            inline const List<label>& superCellForCell() const 
                {return superCellForCell_;}
            
            //- Acces to cellsInSuperCell_
            inline const List<DynamicList<label> >& cellsInSuperCell() const
                {return cellsInSuperCell_;}
       
            //- Number of super cells
            inline label nSuperCells() const 
                {return nSuperCells_;}
            
            //- Number of super cells on current processor
            inline label nSuperCellsOnProcessor() const
                {return nSuperCellsOnProcessor_;}
            
            //- Return the super mesh
            inline const fvMesh& superMesh() const 
                {return superMesh_();}
            
            //- Target number of particles in a super cell
            inline label Npc() const 
                {return Npc_;}          
};
};  // End of namespace Foam

#ifdef NoRepository
    #include "particleNumberControllerTemplates.C"
#endif

#endif


