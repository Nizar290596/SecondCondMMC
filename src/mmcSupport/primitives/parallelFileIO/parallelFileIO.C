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
    Write data to a file using MPI_FILE
    
SourceFiles
    parallelFileIO.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#include "parallelFileIO.H"


void Foam::parallelIO::writeFile
(
    const List<List<scalar>>& data,
    const List<word>& header,
    const fileName& filePath
)
{    
    if (!Pstream::parRun())
    {
        // Revert to writing locally 
        OFstream ofs(filePath);
        writeHeader(ofs, header);
        for (auto& p : data)
        {
            for (auto& e : p)
                ofs << e << token::TAB;
            ofs << nl;
        }
        return;
    }

    // First test implementation collects all information on the master and 
    // writes only from master
    List<List<List<scalar>>> allProcData(Pstream::nProcs());
    allProcData[Pstream::myProcNo()] = data;
    Pstream::gatherList(allProcData);
    
    if (Pstream::master())
    {
        OFstream ofs(filePath);
        writeHeader(ofs, header);
        for (auto& procList : allProcData)
        {
            for (auto& p : procList)
            {
                for (auto& e : p)
                    ofs << e << token::TAB;
                ofs << nl;
            }
        }
    }
}


void Foam::parallelIO::writeFileBinary
(
    const List<List<scalar>>& data,
    const List<word>& header,
    const fileName& filePath
)
{
    // If the fileStatesPtr is not set yet it is now constructed
    if (!fileStatesPtr.valid())
    {
        fileState state;
        // Cast the fileName to the string object it inherits from
        fileStatesPtr.reset
        (
            new std::unordered_map<std::string,fileState> {{filePath,state}}
        );
    }


    auto it = fileStatesPtr->find(filePath);
    if (it == fileStatesPtr->end())
    {
        fileState state;
        fileStatesPtr->insert(std::pair<std::string,fileState>(filePath,state));
        it = fileStatesPtr->find(filePath);
    }

    auto& state = it->second;

    if (Pstream::master() && !state.wroteHeader)
    {
        
        // Check if file already exists
        if (!Foam::exists(filePath))
        {
            OFstream ofs(filePath);
            writeHeader(ofs,header);
            std::ostream& stream = ofs.stdStream();
            state.streamPos = stream.tellp();
        }
        state.wroteHeader = true;
    }

    if (!Pstream::parRun())
    {
        // Revert to writing locally 
        std::ofstream ofs(filePath,std::ios::in|std::ios::out);
        ofs.seekp(state.streamPos);
        
        for (const auto& row : data)
        {
            state.streamPos += row.byteSize();
            ofs.write(reinterpret_cast<const char*>(row.cdata()),row.byteSize());
        }
        
        return;
    }

    // First test implementation collects all information on the master and 
    // writes only from master
    List<List<List<scalar>>> allProcData(Pstream::nProcs());
    allProcData[Pstream::myProcNo()] = data;
    Pstream::gatherList(allProcData);
    
    if (Pstream::master())
    {
        // Revert to writing locally 
        std::ofstream ofs(filePath,std::ios::in|std::ios::out);
        ofs.seekp(state.streamPos);

        for (auto& procList : allProcData)
        {
            for (const auto& row : procList)
            {
                state.streamPos += row.byteSize();
                ofs.write(reinterpret_cast<const char*>(row.cdata()),row.byteSize());
            }
        }
    }
}


void Foam::parallelIO::writeHeader
(
    OFstream& ofs,
    const List<word>& header
)
{
    for (auto& e : header)
        ofs << e << token::TAB;
    ofs << nl;
}


void Foam::parallelIO::setStreamPosBasedOnTime
(
    const double time,
    const label colIndex,
    const fileName& filePath
)
{
    if (Pstream::master())
    {
        // Check that the filePath exists
        if (!Foam::exists(filePath))
            return;

        // Read in the file
        IFstream ifs(filePath);
        DynamicList<List<scalar>> readData;
        
        // Read the header line
        std::string headerLine;
        ifs.getLine(headerLine);
        std::stringstream iss(headerLine);

        // Split the string at seperator
        std::string token;
        label nEntries = 0;
        while(std::getline(iss,token,'\t'))
            nEntries++;


        // Size of one row in bytes
        const label bytesPerRow = nEntries*sizeof(scalar);

        while(ifs.good())
        {
            // Now read binary
            List<scalar> row(nEntries);
            ifs.readRaw(reinterpret_cast<char*>(row.data()),bytesPerRow);
            if (ifs.good())
            {
                readData.append(std::move(row));
            }
        }

        // Loop over the read data and find the correct time index
        forAll(readData,rowI)
        {
            if (readData[rowI][colIndex] >= time)
            {
                // If the fileStatesPtr is not set yet it is now constructed
                if (!fileStatesPtr.valid())
                {
                    fileState state;
                    // Cast the fileName to the string object it inherits from
                    fileStatesPtr.reset
                    (
                        new std::unordered_map<std::string,fileState> {{filePath,state}}
                    );
                }


                auto it = fileStatesPtr->find(filePath);
                if (it == fileStatesPtr->end())
                {
                    fileState state;
                    fileStatesPtr->insert(std::pair<std::string,fileState>(filePath,state));
                    it = fileStatesPtr->find(filePath);
                }

                auto& state = it->second;

                // Set the wroteHeader statement to true
                state.wroteHeader = true;

                // Add the bytes of the header
                state.streamPos = headerLine.size()+1;

                // Set the stream position
                for (label i=0; i < rowI; i++)
                {
                    auto& row = readData[i];
                    state.streamPos += row.byteSize();
                }

                // Break
                break;
            }
        }
    }
}


