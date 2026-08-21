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
    pdfMeshDecomposition

Description
    Arbitrary FV mesh decomposition based on cellDecomposition file of PDF mesh

Procedure
    1) Decompose PDF mesh with any decomposition method with -cellDist option
    to create "cellDecomposition" in /constant/super/ folder
     COMMAND: "decomposePar -region super -cellDist"
    2) Edit system/decomposeParDict file (of FV mesh) and modify the "numbeOfSubdomains"
    to match that of PDF mesh; Modify decomposition method to "manual"
    3) Make this change also:

    manualCoeffs
        {
            dataFile    "cellDecomposition";
        }

    4) run "pdfMeshDecomposition"
    5) Decompose FV mesh ("decomposePar -cellDist"). That's it!!

Author:
    Achinta Varna, UNSW (AV)
\*---------------------------------------------------------------------------*/

#include <string>
#include <sstream>
#include "fvCFD.H"
#include "meshToMesh.H"
#include "bound.H" //added by AV; similar to rhoPimpleFoam;
#include "cell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading cellDecomposition file of PDF mesh from /constant/super folder\n" << endl;
    labelIOList cellDecompositionPDF
    (
        IOobject
        (
            "cellDecomposition",
            runTime.constant(),
            superMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    labelList procIds(mesh.C().size(), -1);

    forAll(mesh.cells(),celli)
    {
        point p = mesh.cellCentres()[celli];

        label superTetFacei = 0;
        label superTetPti = 0;
        label superCelli = 0;

        superMesh.findCellFacePt(p,superCelli,superTetFacei,superTetPti);

        if (superCelli == -1)
        {
            Info << "supercell index: " << superCelli << " for the cell: " << celli << endl;
            superCelli = superMesh.findNearestCell(p);
            Info << "supercell index: " << superCelli << " for the cell: " << celli << endl;
        }

        procIds[celli] = cellDecompositionPDF[superCelli];
    }

    Info<< "Writing cellDecomposition file of FV mesh to /constant folder \n" << endl;
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        procIds
    );

    cellDecomposition.write(); // writing to cellDecomposition file

    Info<< "End pdfMeshDecomposition \n" << endl;

    return 0;
}


// ************************************************************************* //
