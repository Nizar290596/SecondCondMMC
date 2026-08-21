/*
 * radialData.C
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#include "radialData.H"

namespace Foam
{


radialData::radialData(const volVectorField& U_,const fvMesh & mesh_,const Foam::Time &runTime_,const polyPatch& patch,const label& cells)
:
U(U_),
mesh(mesh_),
runTime(runTime_)
{
    fvPatchField<tensor> varPatch = redToPatch(patch, cells);
    redToLine(varPatch);
}

fvPatchField<tensor> radialData::redToPatch(const polyPatch& patch,label cells){        // averages extruded domain to a patch

    volTensorField tmpVar                                                        // initialise temporary Var field for integration
    (
        IOobject
        (
            "tmpVar",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("tmpVar",dimensionSet(0,2,-2,0,0,0,0),tensor::zero),
        "fixedGradient"
    );

    vector Uavg(vector::zero);

    forAll(U, iter){
        tmpVar[iter].xx() = U[iter][0]*U[iter][0];
        tmpVar[iter].xy() = U[iter][0]*U[iter][1];
        tmpVar[iter].xz() = U[iter][0]*U[iter][2];
        tmpVar[iter].yy() = U[iter][1]*U[iter][1];
        tmpVar[iter].yz() = U[iter][1]*U[iter][2];
        tmpVar[iter].zz() = U[iter][2]*U[iter][2];
        Uavg += U[iter];
    }
    Info << "UAvg: " << Uavg/U.size() << endl;

    fvPatchField<tensor> varPatch(tmpVar.boundaryField()[1]);
    fvPatchField<vector> meanPatch(U.boundaryField()[1]);

    forAll(patch, Iter){                                // outer loop  all patch cells
             varPatch[Iter]  = tensor::zero;
             meanPatch[Iter] = vector::zero;
    }


    forAll(patch, Iter){                                // outer loop  all patch cells
        for(int i=0;i<cells;i++){                        // inner loop over all cells along a line
            int pos(Iter+varPatch.size()*i);
             varPatch[Iter]  += tmpVar[pos];
             meanPatch[Iter] += U[pos];
        }
        varPatch[Iter] /= cells;
        meanPatch[Iter] /=cells;
    }

    writeOutVarPatch(varPatch);
    writeOutMeanPatch(meanPatch);

    return varPatch;
}


void radialData::writeOutVarPatch(const fvPatchField<tensor>& patch){

        fileName patchFileTime(runTime.findInstancePath(runTime.findClosestTime(runTime.startTime().value())));
        fileName patchFile(patchFileTime/"inflowUVar");
        OFstream File(patchFile);
        File << "FoamFile\n{\nversion\t2.0;\nformat\tascii;\nclass\tvolTensorField;\nobject\tinflowUvar;\n}\n"; // write stupid header
        File << "dimensions\t[0 2 -2 0 0 0 0];" << endl;
        File << "internalField uniform (0 0 0 0 0 0 0 0 0);" << endl;
        File << "boundaryField\n{\ninletMain\n{\ntype\tfixedValue;\nvalue\tnonuniform List<tensor>\n" << endl;
        File << patch.size() << "\n(";
        forAll(patch, Iter){
            File << patch[Iter] << "\n";
        }
        File << ");\n}\ninletCoflow\n{\ntype\tfixedValue;\nvalue\tuniform (0 0 0 0 0 0 0 0 0);\n}\nfrontAtmosphere\n{\ntype\tfixedValue;\nvalue\tuniform (0 0 0 0 0 0 0 0 0);\n}\nsideAtmosphere\n{\ntype\tzeroGradient;\n";
        File << "}\noutlet\n{\ntype\tzeroGradient;}\n}";

}

void radialData::writeOutMeanPatch(const fvPatchField<vector>& patch){

        fileName patchFileTime(runTime.findInstancePath(runTime.findClosestTime(runTime.startTime().value())));
        fileName patchFile(patchFileTime/"inflowUMean");
        OFstream File(patchFile);
        File << "FoamFile\n{\nversion\t2.0;\nformat\tascii;\nclass\tvolVectorField;\nobject\tinflowUmean;\n}\n"; // write stupid header
        File << "dimensions\t[0 1 -1 0 0 0 0];" << endl;
        File << "internalField uniform (0 0 0);" << endl;
        File << "boundaryField\n{\ninletMain\n{\ntype\tfixedValue;\nvalue\tnonuniform List<vector>\n" << endl;
        File << patch.size() << "\n(";
        forAll(patch, Iter){
            File << patch[Iter] << "\n";
        }
        File << ");\n}\ninletCoflow\n{\ntype\tfixedValue;\nvalue\tuniform (0 0 0);\n}\nfrontAtmosphere\n{\ntype\tfixedValue;\nvalue\tuniform (0 0 0);\n}\nsideAtmosphere\n{\ntype\tzeroGradient;\n";
        File << "}\noutlet\n{\ntype\tzeroGradient;}\n}";

}


void radialData::redToLine(fvPatchField<tensor> varPatch)    // reduces a patch to a line by averaging over phi
{
    // collapse to a single sector
    /*List<scalar> absCor(varPatch.size());
    forAll(varPatch.patch().Cf(),iter){
        absCor[iter] = mag(varPatch.patch().Cf()[iter][0]);
         //mag(varPatch.patch().Cf()[iter][0]);
         //mag(varPatch.patch().Cf()[iter][0]);
    }


    HashSet<scalar>  HashList(absCor.size());
    Foam::List<scalar> myFilteredList();
    label found(0);

    forAll(absCor,Iter){
        if( ! HashList.found(absCor[Iter])){
            found++;
            myFilteredList().resize(found,absCor[Iter]);
            HashList.insert(absCor[Iter]);
        }
    }

    Info << myFilteredList.size() << endl;
    //Info << myFilteredList[0] << endl;*/

    //List<scalar> myFilteredList(HashSet<scalar>(absCor.size()));
    //List<vector> coordinatesCart
    //List<vector> coordinatesRad
    //List<tensor> var

}

}
