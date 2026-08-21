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
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with movingAverage. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "aveSpace/aveSphere/aveSphere.H"

// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * *  //

Foam::aveSphere::aveSphere
(
    const dictionary& dict
)
{
    // Must read length
    radius_ = readScalar(dict.lookup("radius"));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::aveSphere::containsPoint
(
    const point& pCenter,
    const point& p
) const
{
    // Check if point is within LES box by checking the distance 
    const point dist = pCenter-p;

    if (mag(dist) <= radius_)
        return true;
    
    return false;
}

void Foam::aveSphere::printDim() const
{
    // Print dimensions to screen
    Info<< "Averaging sphere " << nl
        << "radius = "   << radius_
        << nl << endl;
}
