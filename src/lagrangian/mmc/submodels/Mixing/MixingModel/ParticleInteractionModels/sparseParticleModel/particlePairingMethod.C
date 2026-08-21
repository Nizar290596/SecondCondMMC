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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#include "particlePairingMethod.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particlePairingMethod::particlePairingMethod(const dictionary& dict)
{
    // If it is a single core run the method is always local
    if (!Pstream::parRun())
    {
        pairingMethod_ = pairingMethod::localPairing;
        return;
    }

    word pairingMethodName = dict.get<word>("pairingMethod");


    if (pairingMethodName == "global")
    {
        pairingMethod_ = pairingMethod::globalPairing;
    }
    else if (pairingMethodName == "local")
    {
        pairingMethod_ = pairingMethod::localPairing;
    }
    else if (pairingMethodName == "subVolumes")
    {
        pairingMethod_ = pairingMethod::subVolumes;
    }
    else
    {
        FatalIOError
            << "Unknown particle pairing method: "<<pairingMethodName<<nl
            << "Possible particle pairing methods are: " << nl
            << token::TAB << "local" << nl
            << token::TAB << "global" <<nl
            << token::TAB << "subVolumes" << exit(FatalIOError); 
    }
}

Foam::Ostream& Foam::operator <<(Ostream& os, const particlePairingMethod& m)
{
    if (m.method() == particlePairingMethod::pairingMethod::localPairing)
        os << "local";
    else if (m.method() == particlePairingMethod::pairingMethod::globalPairing)
        os << "global";
    else if (m.method() == particlePairingMethod::pairingMethod::subVolumes)
        os << "subVolumes";
    
    return os;
}
