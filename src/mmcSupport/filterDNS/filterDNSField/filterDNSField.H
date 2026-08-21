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
    Foam::filterDNSField

Description
    A class providing a filter function to calculate the filtered properties 
    of a field.

    Note: This class requires that the filter mesh is coarser than the DNS mesh.



SourceFiles
    filterDNSField.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef filterDNSField_H
#define filterDNSField_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "meshToMesh.H"
#include "tetOverlapVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                        Class filterDNSField Declaration
\*---------------------------------------------------------------------------*/
class filterDNSField
{
    private:
     
    // Private Member Variables

        //- Reference to the DNS mesh
        const fvMesh& mesh_;

        //- Reference to the filter mesh
        const fvMesh& filterMesh_;

        //- List of cell to filter cell addressing
        //  Each DNS cell is linked to a one or multiple filter cells
        List<DynamicList<label>> cellToFilterCell_;

        //- List of filter cell to cell addressing
        //  Each filter cell is linked to multiple DNS cells
        List<DynamicList<label>> filterCellToCell_;
        
        //- List of cell volume overlap from the DNS cell to the filter cell
        //  A value of 1.0 indicates that the DNS cell is completely in the
        //  filter cell
        List<DynamicList<scalar>> cellToFilterCellVolRatio_;


        //- Processor ID of the filter cell for the given DNS cell ID
        //  List has length of DNS mesh
        //  Can be used to check if the filter cell is on the local processor
        List<DynamicList<label>> cellToFilterCellProcessorID_;

        //- Cells of a local filter cell located on a different processor
        //  For each filter cell a List of Tuple2 is stored which contain 
        //  First: the processor ID of the remote processor
        //  Second: 
        //      Second().first():  the DNS cellID on the remote processor
        //      Second().second(): the volRatio of the DNS cell in the filter cell
        List<DynamicList<Tuple2<label,Tuple2<label,scalar>>>> filterCellToRemoteCells_;

        //- n-cells per filter cell
        //  Access with the cell ID of the filter mesh, has size of the filterMesh
        List<label> nCellsInFilterCell_;


    // Private Member Functions
        
        //- Construct the cellToFilterCell and cellToFilterCellProcessorID list
        void buildCellToFilterLists();

        //- Construct the cellToFilterList for parllel runs
        //  Input: A list of DNS cells that require remote processing stored
        //         as a list of Tuples with the cellID as a label and the 
        //         volume to find on remote processors of that cell as a scalar.
        void buildCellToFilterListsParallel
        (
            const DynamicList<Tuple2<label,scalar>>& DNSCellrequiringRemoteHandling
        );

        //- Set boundary value based on internal cell value
        template<class fieldType>
        void updateBoundary(fieldType& field) const;

        //- Workaround for bug in OpenFOAM to construct the internal boundary
        //  field.
        template<class bFieldType>
        tmp<bFieldType> constructInternalBoundaryField(bFieldType& field) const;

        //- Calculate cell overlap of two cellIDs
        scalar interVol(const label celli, const label filterCellI) const;

        //- Get a list of processor bounding boxes
        List<boundBox> getProcessorBoundingBoxes(const fvMesh& mesh) const;

        //- Calculate the overlap of a bounding box
        scalar overlap
        (
            const label celli,
            const fvMesh& mesh,
            const DynamicList<tetPoints>& tets
        ) const;
    
        //- Calculate the volume overlap of two tets
        void tetTetOverlap
        (
            const tetPoints& tetA,
            const tetPoints& tetB,
            tetPointRef::sumVolOp& insideOp
        ) const;

    public:

        // Constructors

            filterDNSField
            (
                const fvMesh& mesh,
                const fvMesh& filterMesh
            );

            //- Disallow default bitwise copy construct
            filterDNSField(const filterDNSField&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const filterDNSField&) = delete;

        // Filter Functions

            //- Use mean filter 
            template<class fieldType>
            tmp<fieldType> mean(const fieldType& field) const;

            //- Use Favre filter 
            template<class fieldType>
            tmp<fieldType> FavreAverage
            (
                const fieldType& field,
                const volScalarField& rho
            ) const;

            //- Use mean filter to update a field
            //  First DNS field, second the LES field
            template<class fieldType>
            void updateField
            (
                const fieldType& field,
                fieldType& filteredField
            ) const;

            //- Use Favre filter to update a field
            template<class fieldType>
            void updateField
            (
                const fieldType& field,
                const volScalarField& rho,
                fieldType& filteredField
            ) const;

            //- Map filter field to DNS field
            template<class fieldType>
            tmp<fieldType> mapFilterFieldToDNS
            (
                const fieldType& filteredField
            ) const;

            //- Map filter field to DNS field source
            template<class fieldType>
            tmp<fieldType> mapSourceFilterFieldToDNS
            (
                const fieldType& filteredField
            ) const;

        // Const Access

            // Access to the mesh
            const fvMesh& mesh() const {return mesh_;}
            const fvMesh& filterMesh() const {return filterMesh_;}
            
};

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "filterDNSFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
