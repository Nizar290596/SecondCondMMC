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

 
InClass
    Foam::particleTDACChemistryModel
 
Description
    Creates chemistry model instances templated on the type of thermodynamics

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
 
\*---------------------------------------------------------------------------*/
 
#include "makeChemistryModel.H"
 
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
 
#include "particleTDACChemistryModel.H"
#include "thermoPhysicsTypes.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
    // Chemistry moldels based on sensibleEnthalpy
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constGasHThermoPhysics
    //);
 
    makeChemistryModelType
    (
        particleTDACChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );
 
    // Models that can be added if needed

    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constIncompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    incompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    icoPoly8HThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constAdiabaticFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constHThermoPhysics
    //);
 


    ////makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    gasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constIncompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    incompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    icoPoly8HThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constAdiabaticFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constHThermoPhysics
    //);
 
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    constGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    gasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    constIncompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    incompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    icoPoly8HThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    constFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    constAdiabaticFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    psiReactionThermo,
    //    constHThermoPhysics
    //);
 
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    constGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    gasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    constIncompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    incompressibleGasHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    icoPoly8HThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    constFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    constAdiabaticFluidHThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    TDACChemistryModel,
    //    rhoReactionThermo,
    //    constHThermoPhysics
    //);
 
 
    //// Chemistry moldels based on sensibleInternalEnergy
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    gasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constIncompressibleGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    incompressibleGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    icoPoly8EThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constFluidEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constAdiabaticFluidEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    psiReactionThermo,
    //    constEThermoPhysics
    //);
 
 
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    gasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constIncompressibleGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    incompressibleGasEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    icoPoly8EThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constFluidEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constAdiabaticFluidEThermoPhysics
    //);
 
    //makeChemistryModelType
    //(
    //    particleTDACChemistryModel,
    //    rhoReactionThermo,
    //    constEThermoPhysics
    //);
}
