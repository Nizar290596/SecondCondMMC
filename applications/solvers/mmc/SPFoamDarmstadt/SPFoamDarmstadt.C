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
    mmcFoam

Description
    LES / FDF solver for turbulent reacting flows
\*---------------------------------------------------------------------------*/

#include <string>
#include <sstream>
#include "fvCFD.H"
#include "meshToMesh.H"
#include "turbulentFluidThermoModel.H"

#include "basicReactingPopeCloud.H"
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
    #include "createEvalLaplaceMesh.H"
    #include "readGravitationalAcceleration.H"
    Info << g << endl;
    #include "createFields.H"
    Info << "Fields created" << endl;
    #include "createParticles.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createParticleSampling.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting mmcFoam time loop\n" << endl;
    
    scalar Ck_             = 1.5;
    scalar b_vb_           = 1.4;
    scalar gamma_          = 0.5;
    scalar uPrimeCoeff_    = 1.0;
    //scalar fixedsl0_       = 4.09E-01;//SC-Aachen
    //scalar fixedsl0_	     = 0.0421251;//DC-Aachen //Original false
    //scalar fixedsl0_	   = 0.421251; //DC-Aachen adjusted 
    scalar fixedsl0_	   = 3.54E-01; //SC-Darmstadt
    //scalar fixedsl0_	   = 3.66E-01;//DC-Darmstadt
    //scalar fixeddeltal0_   = 3.62E-04;//SC-Aachen
    //scalar fixeddeltal0_   = 0.0020896;//DC-Aachen //Original false
    //scalar fixeddeltal0_   = 0.0004; //DC Aachen adjusted
    scalar fixeddeltal0_   = 3.93E-04; //SC-Darmstadt
    //scalar fixeddeltal0_   = 5.04E-04;//DC-Darmstadt
    scalar eps             = 1e-05;
    
    forAll(dynsl0,cellI)
    {
         dynsl0[cellI]     = fixedsl0_;
         dyndeltal0[cellI] = fixeddeltal0_;
    }

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
	    #include "evaluateLaplacian.H"
	    uPrime = 2.0 * Foam::mag(Foam::pow(LES_delta, 3) * fvc::curl(laplaceSim));
            uPrimek = Foam::sqrt(2.0/3.0 * turbulence->k());

	    Info<< "Calculate lambda function \n" << endl;

	    RL = LES_delta/dyndeltal0;
	    RV = uPrimeCoeff_ * uPrime / dynsl0;
            forAll(vb,cellI) 
            {
                 //scalar a_vb     = 0.6 + 0.2 * Foam::exp(-0.1*RV[cellI]) - 0.2 * Foam::exp(-0.01*RL[cellI]);
                 a_vb[cellI]     = 0.6 + 0.2 * Foam::exp(-0.1*RV[cellI]) - 0.2 * Foam::exp(-0.01*RL[cellI]);
                 //scalar fu       = Foam::max(eps,4.0*Foam::pow(27.0*Ck_/110.0,0.5) * (18.0*Ck_/55.0) * Foam::pow(RV[cellI],2.0));
                 fu[cellI]       = Foam::max(eps,4.0*Foam::pow(27.0*Ck_/110.0,0.5) * (18.0*Ck_/55.0) * Foam::pow(RV[cellI],2.0));
                 //scalar fdelta   = Foam::max(eps,Foam::pow(27.0*Ck_*Foam::pow(3.14,4.0/3.0)/110.0 * (Foam::pow(RL[cellI],4.0/3.0)-1.0),0.5));
                 fdelta[cellI]   = Foam::max(eps,Foam::pow(Foam::max(0.0,27.0*Ck_*Foam::pow(3.14,4.0/3.0)/110.0 * (Foam::pow(RL[cellI],4.0/3.0)-1.0)),0.5));
                 //scalar Resgs    = 4.0 * RL[cellI] * RV[cellI];
                 Resgs[cellI]    = 4.0 * RL[cellI] * RV[cellI];
                 //scalar fre      = Foam::max(eps,Foam::pow(9./55. * Foam::exp(-1.5*Ck_*Foam::pow(3.14,4./3.)/Foam::max(eps,Resgs)),0.5) * Foam::pow(Resgs,0.5));
                 fre[cellI]      = Foam::max(eps,Foam::pow(9./55. * Foam::exp(-1.5*Ck_*Foam::pow(3.14,4./3.)/Foam::max(eps,Resgs[cellI])),0.5) * Foam::pow(Resgs[cellI],0.5));
                 
                 gammafit[cellI] = Foam::pow(Foam::max(0.0,Foam::pow(Foam::max(0.0,Foam::pow(Foam::pow(fu[cellI],-a_vb[cellI])+Foam::pow(fdelta[cellI],-a_vb[cellI]),-1/a_vb[cellI])),-b_vb_) + Foam::pow(fre[cellI],-b_vb_)),-1/b_vb_);
                 vb[cellI]       = Foam::pow(1.0+Foam::min(RL[cellI],gammafit[cellI]*RV[cellI]),gamma_); 
  //               vb[cellI]   = 3.895;a correction factor of 1.4061 must be added to vb when dealing with DC
            }
	    vb.correctBoundaryConditions();
            gammafit.correctBoundaryConditions();

            D = turbulence->nu()/Sc_;
            Dt = turbulence->nut()/Sct_;
            DEff = D + Dt;
//            mu = turbulence->mu();

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

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End mmcFoam " << nl << endl;

    return 0;
}


// ************************************************************************* //
