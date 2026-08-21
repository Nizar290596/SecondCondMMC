/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    CombineBCdata

Description
    Combine inflow data into fewer number of files used by NCIrecycledTimeVaryingMappedFixedValueFvPatchField.

    //skip n time after reading 1 time folder, increase dt of inflow data
    //e.g 0 1 2 3 4 5 6 7 8 9 10 when skipInterval = 2
    //0 3 6 9
    label skipInterval = 0;

    //number of time data per file + 1 end time
    //e.g 0 1 2 3 4 5 6 7 8 9 when nTimes = 3+1
    //{0 1 2 3}{3 4 5 6}{ 6 7 8 9}{9}
    label nTimes = 20+1;

usage:

    Go to target folder e.g. cd xx/xx/xx/constant/boundaryData/inletJet/
    Type command CombineBCdata
    When finished, a new folder named CombinedData containing all BC inflow
    data is created in xx/xx/xx/constant/boundaryData/
    Rename it as needed.
    

Procedure
    1) 

    2) 

    3) 


Author:
    Zhijie Huo (USYD)
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "instantList.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
//    #include "createTime.H"

    
    Info << "Start combining data" << endl;

    mkDir("../CombinedData");

    //find all time folders
    instantList OriginDataTimes = Time::findTimes("./");

    OFstream OriginTimefolders("../CombinedData/OriginTimes");

    OriginTimefolders << OriginDataTimes << endl;

    //skip n time after reading 1 time folder, increase dt of inflow data
    //e.g 0 1 2 3 4 5 6 7 8 9 10 when skipInterval = 2
    //0 3 6 9
    label skipInterval = 0;

    label SkipedTimeSize = std::ceil(
                    scalar(OriginDataTimes.size())/scalar((1+skipInterval))
                                    );

    instantList skipedTimes(SkipedTimeSize);

    for(label i=0;i<SkipedTimeSize;i++)
    {
        skipedTimes[i] = OriginDataTimes[i*(1+skipInterval)];

    }

    OFstream SkipedTimefolders("../CombinedData/SkipTimes");

    SkipedTimefolders << skipedTimes << endl;

    scalar maxdt(0.0);
    label maxI(0);
    for(label i=1;i<skipedTimes.size();i++)
    {
        scalar dt = skipedTimes[i].value()-skipedTimes[i-1].value();
        if(dt > maxdt)
            {
                maxdt = dt;
                maxI = i;
            }

    }

    Info << "max dt = " << maxdt << " at index " << maxI << endl;

    //number of time data per file + 1 end time
    //e.g 0 1 2 3 4 5 6 7 8 9 when nTimes = 3+1
    //{0 1 2 3}{3 4 5 6}{ 6 7 8 9}{9}
    label nTimes = 20+1;

    label CombinedTimeSize = std::ceil(
                    scalar(skipedTimes.size())/scalar((nTimes-1))
                                    );

    Info << "CombinedTimeSize = " << CombinedTimeSize << endl;

    instantList CombinedTimes(CombinedTimeSize);

    List<label> CombinedTimesIndex(CombinedTimeSize+1);

    for(label i=0;i<CombinedTimeSize;i++)
    {
        CombinedTimesIndex[i] = i*(nTimes-1);
        CombinedTimes[i] = skipedTimes[i*(nTimes-1)];

    }

    CombinedTimesIndex.last() = (skipedTimes.size()-1);

    OFstream CombinedTimefolders("../CombinedData/combTimes");

    CombinedTimefolders << CombinedTimes << endl;

    CombinedTimefolders << CombinedTimesIndex << endl;

    forAll(CombinedTimes,i)
    {

        word TimeDir("../CombinedData/"+CombinedTimes[i].name());

        Info << "mkdir " << TimeDir << endl;
        mkDir(TimeDir);

        OFstream combinedU(TimeDir+"/U");

        combinedU << "FoamFile\n{\n"
                    <<"    version     2.0;\n"
                    <<"    format      ascii;\n"
                    <<"    class       dictionary;\n"
                    <<"    object      U;\n"
                    <<"}\n\n\n" << endl;

        for(label indx=CombinedTimesIndex[i];indx<=CombinedTimesIndex[i+1];indx++)
        {
            

            fileName BCdir("./"+skipedTimes[indx].name()+"/U");

            Field<vector> BCdata;

            IFstream(BCdir)() >> BCdata;

            Info << "reading BCdata at " << BCdir << ", size = " << BCdata.size()<< endl;

            combinedU << "\"" << skipedTimes[indx].name() << "\"\n" << BCdata << ";\n\n";

            Info << "written." << endl;

        }




    }
















    Info << "finish combining" << endl;

    return 0;
}


// ************************************************************************* //
