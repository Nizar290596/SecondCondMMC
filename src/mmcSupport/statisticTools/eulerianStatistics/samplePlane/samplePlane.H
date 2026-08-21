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
    Foam::mmcSupport::sampleMesh

Description
    Create a sample plane with function to get the projected face for a given
    point.

    Requires dictionary with:

        rotAxis         (0 0 1);      

        planeType       pointAndNormal;
        triangulate     false;
        pointAndNormalDict
        {   
            point   (0 0 0);
            normal  (0 1.0 0);
        }


SourceFiles
    samplePlane.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef mmcSupport_samplePlane_H
#define mmcSupport_samplePlane_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "kdTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace mmcSupport
{

/*---------------------------------------------------------------------------*\
                        Class samplePlane Declaration
\*---------------------------------------------------------------------------*/
class samplePlane
{
    public:
     
    // Private Member Variables

        const fvMesh& mesh_;
        
        //- Field for coordinates of face centers
        Field<vector2D> coordinates_;
        
        //- Faces
        List<face> faces_;
        
        //- Vertices of the faces
        vectorField points_;
        
        //- Axis for rotating average
        vector rotAxis_;

        //- Normal vector of the plane
        vector normal_;
        
        //- Radial vector of the plane
        //  Perpendicular to normal and axis
        vector rad_;
        
        //- kdTree to quickly look up the closest points 
        autoPtr<kdTree<vector2D>> treePtr_;

        //- Convert point to axial, radial component
        vector2D transform(const vector& p, const bool sign=false) const ;

        //- Combines points and faces of all processors into one plane
        void recombinePointsAndFaces();
    
    public:

        // Constructors

        samplePlane
        (
            const fvMesh& mesh,
            const dictionary& dict
        );
        
        //- Access the size of the sampleMesh as number of cells
        inline label size() const {return coordinates_.size();}
        
        //- Access to the coordinates
        inline const Field<vector2D>& coordinates() const {return coordinates_;}

        //- Return cell index in mesh for given point p
        //  Points are projected onto the sample plane
        label cellIndex(const point& p) const;
        
        //- Write the given field as a vtp file for paraview
        void writeField
        (
            const word& fieldName,
            const Field<scalar>& field
        ) const;

        //- Read plane from file
        //  Return true if file is found and read
        bool read();
            
        //- Write plane to file
        void write() const;     
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mmcSupport
} // End namespace Foam






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
