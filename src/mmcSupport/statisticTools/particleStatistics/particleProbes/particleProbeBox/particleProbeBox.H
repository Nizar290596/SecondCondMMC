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
    Foam::particleparticleProbeBox

Description
    Probing volume representing a box constructed from two points
    
Usage
    
    dictionary
    {
        type particleProbeBox;
        p1      (1 0 0);
        p2      (2 0 0);
    } 
    

SourceFiles
    particleparticleProbeBox.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020>

\*---------------------------------------------------------------------------*/

#ifndef particleProbeBox_H 
#define particleProbeBox_H



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "particleProbes.H"
#include "volFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class particleProbeBox Declaration
\*---------------------------------------------------------------------------*/

class particleProbeBox
:
    public particleProbes
{
    private:
        
        point p1_;
        point p2_;

public:

    //- Runtime type information
    TypeName("particleProbeBox");

    // Constructor

        //- Construct from dictionary 
        particleProbeBox(const dictionary& dict);
        
    // Public Member Functions

        //- Return true if the particle is within the probing space
        virtual bool containsParticle(const particle& p) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
