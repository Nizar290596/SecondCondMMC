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
    Foam::ParticleInCell
 
Description
    Couples thermophysical properties of the particles with the LES field by
    averaging all particles within the super cell of the current cell processed
    

SourceFiles
    ParticleInCell.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef ParticleInCell_H
#define ParticleInCell_H

#include "ThermoPhysicalCouplingModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ParticleInCell Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleInCell
:
    public ThermoPhysicalCouplingModel<CloudType>
{
    // Private data
    const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("ParticleInCell");

    //- Type of particle the cloud was instantiated for
    typedef typename CloudType::particleType particleType;
    
    // Constructors

        //- Construct from dictionary
        ParticleInCell
        (
            const dictionary& dict, 
            CloudType& owner, 
            const mmcVarSet& Xi
        );

        //- Destructor
        virtual ~ParticleInCell();

    // Sources

        //- Equivalent temperature and species target values
        virtual void EqvETargetValues
        (
            const mmcVarSet& Xi,
            PtrList<volScalarField>& YEqvETarget,
            volScalarField& TEqvETarget
        );

  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleInCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
