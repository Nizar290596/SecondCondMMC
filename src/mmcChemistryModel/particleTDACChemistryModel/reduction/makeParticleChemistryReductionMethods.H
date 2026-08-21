/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#ifndef makeParticleChemistryReductionMethods_H
#define makeParticleChemistryReductionMethods_H

#include "particleChemistryReductionMethod.H"

#include "noParticleChemistryReduction.H"
#include "pDAC.H"
#include "pDRG.H"
#include "pDRGEP.H"
#include "pEFA.H"
#include "pPFA.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 
#define makeParticleChemistryReductionMethod(SS, Comp, Thermo)                 \
                                                                               \
    typedef particleChemistryReductionMethods::SS<Comp, Thermo>                \
        SS##Comp##Thermo;                                                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    particleChemistryReductionMethod<Comp, Thermo>::                           \
        adddictionaryConstructorToTable<SS##Comp##Thermo>                      \
        add##particleChemistryReductionMethods##SS##Comp##Thermo##ConstructorToTable_;
 
 
#define makeParticleChemistryReductionMethods(CompChemModel, Thermo)           \
                                                                               \
    typedef particleChemistryReductionMethod<CompChemModel, Thermo>            \
        particleChemistryReductionMethod##CompChemModel##Thermo;               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        particleChemistryReductionMethod##CompChemModel##Thermo,                       \
        (word(particleChemistryReductionMethod##CompChemModel##Thermo::typeName_()) +  \
        '<' + word(CompChemModel::typeName_()) + "," + Thermo::typeName() + '>'\
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        particleChemistryReductionMethod##CompChemModel##Thermo,               \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeParticleChemistryReductionMethod                                       \
    (                                                                          \
        pNone,                                                                 \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
 makeParticleChemistryReductionMethod                                          \
    (                                                                          \
        pDAC,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistryReductionMethod                                       \
    (                                                                          \
        pDRG,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistryReductionMethod                                       \
    (                                                                          \
        pDRGEP,                                                                \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistryReductionMethod                                       \
    (                                                                          \
        pEFA,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeParticleChemistryReductionMethod                                       \
    (                                                                          \
        pPFA,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
