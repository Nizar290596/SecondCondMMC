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
    AxisymAv

Description
    Circumferential average

Author:
    Zhijie Huo (USYD)
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "instantList.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "vector.H"
//#include <algorithm>

#include <string>
#include <sstream>
#include "fvCFD.H"

#include "fvOptions.H"
#include "cell.H"

using namespace Foam;


label findIntervalIndex( const scalar& val, const List<scalar>& pointList)
{

    label left = 0;
    label right = pointList.size()-1;

    while(left < right)
    {
        label midpoint = std::floor( left + (right - left) / 2 );

        if(pointList[midpoint] < val)
        {
            left = midpoint+1;
        }
        else
        {
            right = midpoint;
        }

    }

    return left-1;

}



int main(int argc, char *argv[])
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // execution
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    //#include "createFvOptions.H"

    //#include "createControl.H"
    //#include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField meshCRadius
    (
        IOobject
        (
            "meshCRadius",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.C().component(vector::X)*0
    );


    forAll(meshCRadius,celli)
    {
        scalar Rtemp2 = pow(mesh.C()[celli].component(vector::X),2) 
                            + pow(mesh.C()[celli].component(vector::Y),2);

        meshCRadius[celli] = std::sqrt(Rtemp2);
    }

    //Group cells by Z coordinate.

    Info << "Grouping Z coordinates" << endl;

    List<scalar> AllZgridPonits(mesh.points().component(vector::Z));

    std::sort(AllZgridPonits.begin(),AllZgridPonits.end());

    List<scalar>::iterator itZ;

    itZ = std::unique(AllZgridPonits.begin(),AllZgridPonits.end());

    AllZgridPonits.resize(std::distance(AllZgridPonits.begin(),itZ));

   // Info << "AllZgridPonits = " << AllZgridPonits << endl;

    //container of grouped cell centre ID, size=points-1
    List< List<label> > Layers_CellID(AllZgridPonits.size()-1, List<label>());

    forAll(mesh.C(),celli)
    {
        label Layer_i = findIntervalIndex(mesh.C()[celli].component(vector::Z), AllZgridPonits);

        Layers_CellID[Layer_i].append(celli);

    }

/*
        Info << "mesh.C()[Layers_CellID[1][5]] = " << mesh.C()[Layers_CellID[1][5]] << endl;
        Info << "mesh.C()[Layers_CellID[1][15]] = " << mesh.C()[Layers_CellID[1][15]] << endl;
        Info << "mesh.C()[Layers_CellID[1][25]] = " << mesh.C()[Layers_CellID[1][25]] << endl;
        Info << "mesh.C()[Layers_CellID[1][35]] = " << mesh.C()[Layers_CellID[1][35]] << endl;
        Info << "mesh.C()[Layers_CellID[1][45]] = " << mesh.C()[Layers_CellID[1][45]] << endl;

        Info << "mesh.C()[Layers_CellID[0][50]] = " << mesh.C()[Layers_CellID[0][50]] << endl;

        Info << "mesh.C()[Layers_CellID.last()[5]] = " << mesh.C()[Layers_CellID.last()[5]] << endl;
        Info << "mesh.C()[Layers_CellID.last()[15]] = " << mesh.C()[Layers_CellID.last()[15]] << endl;
        Info << "mesh.C()[Layers_CellID.last()[25]] = " << mesh.C()[Layers_CellID.last()[25]] << endl;
*/


    //group cellID in each layer by radius
    Info << "Grouping R coordinates" << nl<< endl;
    
    List<scalar> AllRgridPonits(mesh.points().component(vector::X));

    //Info << "AllRgridPonits.size() = " << AllRgridPonits.size() << endl;

    scalar Radius_max = max(AllRgridPonits);

    Info << "Radius_max = " << Radius_max << endl;

    scalar res_R = 0.0001; //m

    label divPointsNum = std::ceil(Radius_max / res_R) +1;

    //radius grid point for goruping
    List<scalar> divGridPoints(divPointsNum);

    for(label Di = 0; Di < divPointsNum; Di++)
    {
        divGridPoints[Di] = res_R * Di;

    }

    //Info << "divGridPoints[0] = " << divGridPoints[0] << endl;
    //Info << "divGridPoints.last() = " << divGridPoints.last() << endl;


    label NumLayer = Layers_CellID.size();

    label NumRadDiv = std::ceil(Radius_max / res_R);
    
    List<List< List<label> >>   GroupedCellID(NumLayer, List< List<label> >(NumRadDiv));
    //List<List< List<label> >>   GroupedCellID(1, List< List<label> >(NumRadDiv));

    forAll(GroupedCellID, layerI)
    {

        List<label> ith_layer(Layers_CellID[layerI]);

        List< List<label>>  ith_layerTemp(NumRadDiv);

        forAll(ith_layer, Radi)
        {
            label RadGroupi = findIntervalIndex( meshCRadius[ith_layer[Radi]], divGridPoints); 

            //Info << "meshCRadius[ith_layer[Radi]] = " << meshCRadius[ith_layer[Radi]] << endl;
            //Info << "RadGroupi = " << RadGroupi << endl;


            ith_layerTemp[RadGroupi].append(ith_layer[Radi]);

        }


        //Info << "ith_layerTemp[0] = " << ith_layerTemp[0] << endl;

        GroupedCellID[layerI]=(ith_layerTemp);

        //Info << "GroupedCellID[layerI] = " << GroupedCellID[layerI] << endl;


    }

/*
    Info << "GroupedCellID.size() = " << GroupedCellID.size() << endl;

    Info << "GroupedCellID[0].size() = " << GroupedCellID[0].size() << endl;

    //Info << "GroupedCellID[0][0].size() = " << GroupedCellID[0][0].size() << endl;

    //Info << "GroupedCellID[0] = " << GroupedCellID[0] << endl;


    forAll(GroupedCellID[0][0],ii)
    {Info << "meshCRadius[GroupedCellID[0][0]["<<ii<<"]] = " << meshCRadius[GroupedCellID[0][0][ii]] << endl;}


    forAll(GroupedCellID[457][0],ii)
    Info << "meshCRadius[GroupedCellID[457][0]["<<ii<<"]] = " << meshCRadius[GroupedCellID[457][0][ii]] << endl;
*/





    forAll(timeDirs, timeI)
    {
       	runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        forAll(AvFields,Yi)
        {

            forAll(GroupedCellID,layerI)
            {
                forAll(GroupedCellID[layerI], ri)
                {

                    if(GroupedCellID[layerI][ri].size())
                    {   
                        scalar AvAccum = 0;
                        label  AvDataNum = 0;

                        forAll(GroupedCellID[layerI][ri],ci)
                        {
                            label cellIDi = GroupedCellID[layerI][ri][ci];
                
                            AvAccum += AvFields[Yi].internalField()[cellIDi];

                            AvDataNum ++;

                        }
                        
                        scalar IntervalAv = AvAccum/AvDataNum;

                        forAll(GroupedCellID[layerI][ri],ci)
                        {
                            label cellIDi = GroupedCellID[layerI][ri][ci];
                
                            CircAvFields[Yi][cellIDi] = IntervalAv;

                        }


                    }



                }


            }

            CircAvFields[Yi].write();

            //Info << "Tav cell number = " << TAv.internalField().size() << endl;

            //Info << "T_max = " << max(TAv) << endl;

            //Info << "mesh.C cell number = " << mesh.C().size() << endl;


        }
    }

}

























