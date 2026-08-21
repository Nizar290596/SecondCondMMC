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

#include "mmcStatusMessage.H"


Foam::word Foam::mmcStatusMessage::gitBranch_=GIT_BRANCH;
Foam::word Foam::mmcStatusMessage::gitHash_=GIT_BUILD;
Foam::fileName Foam::mmcStatusMessage::mmcBasePath_=MMC_PROJECT_USER_DIR;



Foam::fileName Foam::mmcStatusMessage::solverBasePath(const fileName& pathToSolver) 
{
    return mmcBasePath_/pathToSolver;
}


void Foam::mmcStatusMessage::print() 
{
    const word build(gitBranch_ + " " + gitHash_);
    
    const label strSize = 62-build.length(); 
   
    fileName installPath = mmcBasePath_;

    // Check that the length of the file path is not too long
    if (mmcBasePath_.length() > 59)
    {
        installPath = mmcBasePath_.substr(0,59); 
    }
    

    Info << nl
         << "/*---------------------------------------------------------------------------*\\"<< nl
         << "                                       8888888888                              " << nl 
         << "                                       888                                     " << nl 
         << "                                       888                                     " << nl 
         << "  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  " << nl 
         << "  888 \"888 \"88b 888 \"888 \"88b d88P\"    888     d88\"\"88b     \"88b 888 \"888 \"88b " << nl 
         << "  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 " << nl 
         << "  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 " << nl 
         << "  888  888  888 888  888  888  \"Y8888P 888      \"Y88P\"  \"Y888888 888  888  888 " << nl 
         << "-------------------------------------------------------------------------------"<<nl
         << "|   Build:       "  << build <<  setw(strSize)                               << "|" << nl
         << "|   Installed in "  << installPath << setw(60-installPath.length())          << "|" << nl 
         << "\\*---------------------------------------------------------------------------*/" << nl
         << endl;
}


