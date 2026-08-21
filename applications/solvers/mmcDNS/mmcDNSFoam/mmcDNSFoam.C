/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mmcDNSFoam

Group
    grpCombustionSolvers

Description
    Solves the reactive scalars on the Eulerian grid based on rhoReactiveFoam.
    The mmc clouds receive information and are passively tracked with the
    simulation.

\*---------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "filterDNSField.H"
#include "filteredTurbulence.H"

#include "basicReactingPopeCloud.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cell.H"

//MMC Item
#include "mmcVarSet.H"
#include "multivariateScheme.H"
#include "mmcStatusMessage.H"
#include "particleStatistics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver Carrier-Phase DNS as Reference for MMCLES beased on rhoReactingFoam"
        " thermodynamics package."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    // Print a status message
    mmcStatusMessage::print();
    #include "createTimes.H"
    #include "createMesh.H"
    #include "createFilterMesh.H"

    #include "readGravitationalAcceleration.H"
    Info << g << endl;
    
    #include "createControl.H"
    #include "initContinuityErrs.H"
    #include "readPDFDictionarySettings.H"

    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFilteredFields.H"
    #include "updateFilteredFields.H"

    #include "createParticles.H"
    #include "createRhoUfIfPresent.H"
    #include "createTimeControls.H"

    #include "createCloudSampling.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        // Set the filterTime step to DNS time step -- for consistency
        filterTime.setDeltaT(runTime.deltaT());

        // Update DNS and filter time
        runTime++; // We need the time 0 to inject ther particle
        filterTime++;


        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H" // we take the equation from MMCLES

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            
            // Update the diffusivity
            D = turbulence->nu()/Sc_;

            #include "XiEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        // Update all filtered fields
        #include "updateFilteredFields.H"

        #include "moveParticles.H"    

        runTime.write();
        filterTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;   
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
