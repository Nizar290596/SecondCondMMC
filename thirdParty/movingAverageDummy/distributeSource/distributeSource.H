/*---------------------------------------------------------------------------*\
    888b     d888                                  d8888                   
    8888b   d8888                                 d88888                   
    88888b.d88888                                d88P888                   
    888Y88888P888  .d88b.  888  888  .d88b.     d88P 888 888  888  .d88b.  
    888 Y888P 888 d88""88b 888  888 d8P  Y8b   d88P  888 888  888 d8P  Y8b 
    888  Y8P  888 888  888 Y88  88P 88888888  d88P   888 Y88  88P 88888888 
    888   "   888 Y88..88P  Y8bd8P  Y8b.     d8888888888  Y8bd8P  Y8b.     
    888       888  "Y88P"    Y88P    "Y8888 d88P     888   Y88P    "Y8888  
-------------------------------------------------------------------------------                                                                                                                                                    
License
    This file is part of movingAverage.

    movingAverage is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    movingAverage is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with movingAverage.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::distributeSource

Description
    This is a dummy class if the movingAverage library cannot be loaded

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024

\*---------------------------------------------------------------------------*/

#ifndef distributeSource_H 
#define distributeSource_H

#include "aveSpace/aveSpace/aveSpace.H"
#include "fvMesh.H"
#include "regIOobject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Dummy class for movingLESCellManager
\*---------------------------------------------------------------------------*/
class movingLESCellManager
{
    public:
        // Dummy default constructor
        movingLESCellManager() = default;
};

/*---------------------------------------------------------------------------*\
                  Class distributeSource Declaration
\*---------------------------------------------------------------------------*/
class distributeSource
:
    public regIOobject
{
    private:

        movingLESCellManager manager_;
    
    public:

    // Constructor
       
        distributeSource(const fvMesh& mesh, const aveSpace& box)
        : regIOobject
        (
            IOobject
            (
                "distributeSourceDummy",
                mesh.time().constantPath(),
                mesh
            )
        ) {}

        distributeSource(const fvMesh& mesh, const movingLESCellManager& manager)
        : regIOobject
        (
            IOobject
            (
                "distributeSourceDummy",
                mesh.time().constantPath(),
                mesh
            )
        ) {}


        //- Disallow default bitwise copy construct
        distributeSource(const movingLESCellManager&) = delete;

        //- Disallow default bitwise assignment
        distributeSource& operator=(const movingLESCellManager&) = delete;

  // Member Functions

    // Access

        //- Floating LES cell manager
        inline const movingLESCellManager& manager() const
        {
            return manager_;
        }

    // Calculations

        //- Distribute the cell values to the given cells in the defined 
        //- LES box returns again a volField with the distributed values
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> distribute
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            FatalError << "Calling the dummy class -- not implemented"
                << exit(FatalError);
            tmp<GeometricField<Type, fvPatchField, volMesh>> empty;
            return empty;
        }
        
        //- Distribute the cell values to the given cells in the defined 
        //- LES box returns again a volField with the distributed values
        template<class Type>
        tmp<DimensionedField<Type, volMesh>> distribute
        (
            const DimensionedField<Type, volMesh>& vf
        ) const
        {
            FatalError << "Calling the dummy class -- not implemented"
                << exit(FatalError);
            tmp<DimensionedField<Type, volMesh>> empty;
            return empty;   
        }

    // Write Data
        
        //- Write data inherited from regIOobject
        bool writeData(Ostream& os) const {return true;}
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
