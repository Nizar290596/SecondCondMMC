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

Description
    Create the different chmistry ode solver types 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/
 
#include "makeParticleChemistrySolverTypes.H"
 
#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
    // To add additional solver types all make* files have to be adjusted
    // E.g. for a constGasHThermoPhysics case the 
    // * makeParticleChemistryModels.C
    // * makeParticleTDACChemistryModels.C
    // * makeParticleChemistryReductionMethods.C
    // * makeParticleChemistryTabulationMethods.C
    // have to be adjusted as well with a model for constGasHThermoPhysics
    // otherwise a undefined reference error occurs


    // Chemistry solvers based on sensibleEnthalpy
    makeParticleChemistrySolverTypes
    (
        psiReactionThermo, 
        gasHThermoPhysics
    );
}
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


