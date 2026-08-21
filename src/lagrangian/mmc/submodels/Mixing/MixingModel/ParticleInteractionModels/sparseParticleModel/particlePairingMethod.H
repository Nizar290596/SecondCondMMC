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
    Foam::particlePairingMethod
 
Description
    Provides a more descriptive user interface to select different particle
    pariring methods. The actual particle pairing is done in mixParticleModel

SourceFiles
    particlePairingMethod.C
    mixParticleModel.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef particlePairingMethod_H
#define particlePairingMethod_H

#include "dictionary.H"
#include "Switch.H"
#include "Pstream.H"
#include "Ostream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declaration
class particlePairingMethod;

// Write information of particle model to stream
Ostream& operator << (Ostream& os, const particlePairingMethod& m);


/*---------------------------------------------------------------------------*\
                Class particlePairingMethod Declaration
\*---------------------------------------------------------------------------*/

class particlePairingMethod
{
    public: 

        enum pairingMethod
        {
            localPairing,
            subVolumes,
            globalPairing
        };

    private:

        pairingMethod pairingMethod_;


    public:

        // Constructor

            //- Construct from dictionary
            particlePairingMethod(const dictionary& dict);

        // Access

            //- Return selected particle pairing method
            inline pairingMethod method() const 
                {return pairingMethod_;}

            //- Is it a global method
            inline Switch global() const 
                {return pairingMethod_ == pairingMethod::globalPairing;}

};
} // End of namespace Foam

#endif

// ************************************************************************* //