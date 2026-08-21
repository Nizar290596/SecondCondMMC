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

#include "particleProbeBox.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleProbeBox, 0);
    addToRunTimeSelectionTable(particleProbes, particleProbeBox, particleProbeType);
}

// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * *  //

Foam::particleProbeBox::particleProbeBox
(
    const dictionary& dict
)
:
    particleProbes(dict)
{
    // Read the two points
    point p1Temp = dict.get<point>("p1");
    point p2Temp = dict.get<point>("p2");
    
    for (label i=0; i < 3; i++)
    {
        if (p1Temp[i] > p2Temp[i])
        {
            p1_[i] = p2Temp[i];
            p2_[i] = p1Temp[i];
        }
        else
        {
            p1_[i] = p1Temp[i];
            p2_[i] = p2Temp[i];
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::particleProbeBox::containsParticle
(
    const particle& p
) const
{
    point pos = p.position();
    
    if
    (
        pos.x() >= p1_.x() && pos.x() <= p2_.x()
     && pos.y() >= p1_.y() && pos.y() <= p2_.y()
     && pos.z() >= p1_.z() && pos.z() <= p2_.z()
    )
        return true;
    
    return false;
}
