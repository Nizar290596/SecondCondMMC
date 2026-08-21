/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
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
    mmcPremixedFoam

Description
    LES / FDF solver for turbulent reacting flows
\*---------------------------------------------------------------------------*/

#include <string>
#include <sstream>
#include "fvCFD.H"
#include "meshToMesh.H"
#include "turbulentFluidThermoModel.H"

#include "basicPremixedReactingPopeCloud.H"
#include "psiReactionThermo.H"
#include "LESdelta.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "cell.H"

#include "mmcVarSet.H"
#include "multivariateScheme.H"
#include "mmcStatusMessage.H"
#include "particleStatistics.H"

#include "interpolateFlameletTable.H" 
#include "efficiencyFunction.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    // Print a status message
    mmcStatusMessage::print();
    #include "createTime.H"
    #include "createMesh.H"
    #include "createEvalLaplaceMesh.H" 
    #include "readGravitationalAcceleration.H"
    Info << g << endl;
    #include "createFields.H"
    #include "readFlameletTable.H" 
    #include "partInitRelease_hYEqvE_Z42.H" 
    #include "createParticlesPremixed.H"
    #include "flamelet_hYEqvE_Z42.H" 
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createParticleSampling.H"
    #include "readMixTimeTable.H" 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting mmcPremixedFoam time loop" << nl << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        // --- PIMPLE loop
        while (pimple.loop())
        {
            #include "UEqnTurbForcing.H"
            
            D = turbulence->nu()/Sc_;
            Dt = turbulence->nut()/Sct_;
            DEff = D + Dt;
            mu = turbulence->mu();
            
            #include "flameletClosure_pvSource.H" 
            #include "calcATFproperties.H" 
            #include "XiEqnATF.H"
            #include "flamelet_hYEqvE_Z42.H" 

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        gradRho = fvc::grad(rho*DEff);

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;

        #include "moveParticlesPremixed.H"

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End mmcPremixedFoam " << nl << endl;

    return 0;
}


// ************************************************************************* //
