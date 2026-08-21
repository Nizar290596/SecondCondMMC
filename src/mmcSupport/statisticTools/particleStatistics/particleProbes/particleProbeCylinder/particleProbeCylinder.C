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

#include "particleProbeCylinder.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleProbeCylinder, 0);
    addToRunTimeSelectionTable(particleProbes, particleProbeCylinder, particleProbeType);
}

// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * *  //

Foam::particleProbeCylinder::particleProbeCylinder
(
    const dictionary& dict
)
:
    particleProbes(dict),
    p1_(dict.get<point>("p1")),
    p2_(dict.get<point>("p2")),
    radius_(dict.get<scalar>("radius"))
{
    // Calculate the normal direction
    axis_ = (p2_- p1_);
    axis_ /= mag(axis_);

    length_ = mag(p2_-p1_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::particleProbeCylinder::containsParticle
(
    const particle& p
) const
{
    point pos = p.position();
    
    // Calculate if it is in the axial range of cylinder
    scalar axialLength = (pos - p1_) & axis_;

    if (axialLength > 0 && axialLength < length_)
    {
        // Check if it is in radial dimension 
        scalar rad = mag(axis_^(pos-p1_));
        if (rad < radius_)
            return true;
    }
    
    return false;
}
