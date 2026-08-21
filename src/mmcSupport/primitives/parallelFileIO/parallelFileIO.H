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
 
Description
    Write data to a file taking care of parallel runs 
    
SourceFiles
    parallelFileIO.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#ifndef parallelFileIO_H
#define parallelFileIO_H

#include "List.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Pstream.H"
#include "IOdictionary.H"
#include <unordered_map>

namespace Foam
{
namespace parallelIO
{
    struct fileState
    {
        std::streamsize streamPos = 0;
        bool wroteHeader = false;
    };

    static autoPtr<std::unordered_map<std::string,fileState>> fileStatesPtr;

    void writeFile
    (
        const List<List<scalar>>& data,
        const List<word>& header,
        const fileName& filePath
    );

    void writeFileBinary
    (
        const List<List<scalar>>& data,
        const List<word>& header,
        const fileName& filePath
    );

    
    void writeHeader
    (
        OFstream& ofs,
        const List<word>& header
    );

    // Set the streamPos index based on the given time value
    // This is required for a restart from a previous time point, which is 
    // not the last time, and collated fileIO.
    // E.g., a simulation runs from 0 to 0.5s and then is restarted at 0.4s.
    // This avoids duplicate entries for the time between 0.4s and 0.5s.
    void setStreamPosBasedOnTime
    (
        const double time,
        const label colIndex,
        const fileName& filePath
    );
}

};  // End of namespace Foam

#endif
