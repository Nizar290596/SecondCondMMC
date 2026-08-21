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

Class
    Foam::particleChemistrySolverTypes
 
Description
    Macro to create the different chmistry ode solver types 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/
 
#ifndef makeChemistrySolverTypes_H
#define makeChemistrySolverTypes_H
 
#include "chemistrySolver.H"
 
#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"
 
#include "particleChemistryModel.H"
#include "particleTDACChemistryModel.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#define makeParticleChemistrySolverType(SS, Comp, Thermo)                      \
                                                                               \
    typedef SS<particleChemistryModel<Comp, Thermo>> SS##Comp##Thermo;         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(particleChemistryModel<Comp, Thermo>::typeName_()) + "<"\
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_;                                                                      \
                                                                               \
    typedef SS<particleTDACChemistryModel<Comp, Thermo>>                       \
        TDAC##SS##Comp##Thermo;                                                \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        TDAC##SS##Comp##Thermo,                                                \
        (#SS"<" + word(particleTDACChemistryModel<Comp, Thermo>::typeName_())  \
        + "<" + word(Comp::typeName_()) + "," + Thermo::typeName() +           \
        ">>").c_str(), 0                                                       \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<TDAC##SS##Comp##Thermo>                \
        add##TDAC##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel\
##Comp##Table_;
 
 
#define makeParticleChemistrySolverTypes(Comp, Thermo)                         \
                                                                               \
    makeParticleChemistrySolverType                                            \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistrySolverType                                            \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistrySolverType                                            \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //
