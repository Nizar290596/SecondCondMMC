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
    mmcSprayFoam

    Includes the DropletSpray cloud and particles, for spray combustion.

Description
    LES / FDF solver for turbulent reacting flows
\*---------------------------------------------------------------------------*/

#include <string>
#include <sstream>
#include "fvCFD.H"
#include "meshToMesh.H"
#include "turbulentFluidThermoModel.H"

#include "basicReactingPopeCloud.H"
#include "basicDropletSprayThermoCloud.H"
#include "psiReactionThermo.H"
#include "LESdelta.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "cell.H"

#include "mmcVarSet.H"
#include "multivariateScheme.H"
#include "mmcStatusMessage.H"
#include "particleStatistics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    // Print a status message
    mmcStatusMessage::print();
    #include "createTime.H"
    #include "createMesh.H"

    #include "readGravitationalAcceleration.H"
    Info << g << endl;
    #include "createFields.H"
    #include "createParticles.H"
    #include "createFuelCloud.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createParticleSampling.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting mmcFoam time loop" << nl << endl;

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
            #include "UEqn.H"

            D = turbulence->nu()/Sc_;
            Dt = turbulence->nut()/Sct_;
            DEff = D + Dt;
            mu = turbulence->mu();

            #include "XiEqn.H"
            #include "hYEqvE_Eqn.H"

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

        #include "moveParticles.H"

        #include "moveFuelParticles.H"

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End mmcFoam " << nl << endl;

    return 0;
}


// ************************************************************************* //
