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
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Application
    writeSubVolumes

    Write the generated sub volume lists as processor ids to file

Description
    LES / FDF solver for turbulent reacting flows
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mixingSubVolumes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    // Only run in parallel 
    if (!Pstream::parRun())
        FatalError << "Test must be run in parallel" << exit(FatalError);

    const IOdictionary cloudPropertiesDict
    (
        IOobject
        (
            "cloudProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    dictionary subModelDict = cloudPropertiesDict.subDict("subModels");

    dictionary mixingDict = subModelDict.subDict("MMCcurlCoeffs");

    const scalar ri = mixingDict.get<scalar>("ri");
    const scalar Xii = mixingDict.subDict("Ximi").get<scalar>("fm");

    // Test the constructor
    mixingSubVolumes mixSubVolumes(mesh,mixingDict,ri,Xii);

    // Write subVolumeSets to visualize in paraview
    mixSubVolumes.writeSubVolumeSets();
}

// ************************************************************************* //
