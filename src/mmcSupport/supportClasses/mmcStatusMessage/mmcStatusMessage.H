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
    Foam::mmcStatusMessage
 
Description
    Prints a welcome message to the screen and provides access 
    to the solver base path required for the mmcVarSet class

SourceFiles
    mmcStatusMessage.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/




#ifndef mmcStatusMessage_H
#define mmcStatusMessage_H

#ifndef GIT_BRANCH
    #define GIT_BRANCH "NaN"
#endif

#ifndef GIT_BUILD
    #define GIT_BUILD "NaN"
#endif

#ifndef MMC_PROJECT_USER_DIR
    #define  MMC_PROJECT_USER_DIR "Solver_path_not_set_correctly"
#endif


#include "IOobject.H"
#include "IOmanip.H"
#include <string>

namespace Foam
{

class mmcStatusMessage
{
    private:
    
        static word gitBranch_;
        static word gitHash_;
       
        static fileName mmcBasePath_;
        
    public:

        //- Constructor 
        mmcStatusMessage() = default;
        
        //- Returns the path to the solver path
        static fileName solverBasePath(const fileName& pathToSolver);

        //- Prints a welcome message with mmcFoam banner 
        //- to the screen
        static void print();
};

};  // End of namespace Foam

#endif
