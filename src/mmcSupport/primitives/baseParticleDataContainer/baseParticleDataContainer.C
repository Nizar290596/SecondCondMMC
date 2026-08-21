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

#include "baseParticleDataContainer.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * * 
Foam::baseParticleDataContainer::baseParticleDataContainer(Istream& is)
{
    read(is);
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * *

void Foam::baseParticleDataContainer::append(baseParticleDataContainer& container)
{
    labelValues_.append(container.labelValues_);
    scalarValues_.append(container.scalarValues_);
    scalarFields_.append(container.scalarFields_);
    
    // the position is not modified
}


// * * * * * * * * * * * * * * * IO Functions * * * * * * * * * * * * * * * * *

Foam::Istream& Foam::baseParticleDataContainer::read
(
    Istream& is
)
{
    labelValues_.clear();
    is >> labelValues_;
    
    scalarValues_.clear();
    is >> scalarValues_;
    
    scalarFields_.clear();
    is >> scalarFields_;
    
    is >> pos_;

    is >> local_;

    return is;
}

Foam::Ostream& Foam::baseParticleDataContainer::write
(
    Ostream& os
) const
{
    // When the particle container is transferred or written it is no longer
    // local
    bool localParticle = false;

    os << labelValues_ << token::SPACE
       << scalarValues_ << token::SPACE
       << scalarFields_ << token::SPACE
       << pos_ 
       << localParticle << endl;
    return os;
}


Foam::Istream& Foam::operator >>(Istream& is, baseParticleDataContainer& eField)
{
    return eField.read(is);
}


Foam::Ostream& Foam::operator <<(Ostream& os, const baseParticleDataContainer& eField)
{
    return eField.write(os);
}
