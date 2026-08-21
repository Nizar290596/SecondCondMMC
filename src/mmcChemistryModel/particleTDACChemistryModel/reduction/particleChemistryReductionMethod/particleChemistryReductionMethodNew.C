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

#include "particleChemistryReductionMethod.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::particleChemistryReductionMethod<CompType, ThermoType>>
Foam::particleChemistryReductionMethod<CompType, ThermoType>::New
(
    const IOdictionary& dict,
    particleTDACChemistryModel<CompType, ThermoType>& chemistry
)
{
     const dictionary& reductionDict = dict.subDict("reduction");
  
     const word methodName(reductionDict.get<word>("method"));
  
     Info<< "Selecting chemistry reduction method " << methodName << endl;
  
     const word methodTypeName
     (
         methodName
       + '<' + CompType::typeName + ',' + ThermoType::typeName() + '>'
     );
  
     auto cstrIter = dictionaryConstructorTablePtr_->cfind(methodTypeName);
  
     if (!cstrIter.found())
     {
         constexpr const int nCmpt = 7;
  
         wordList thisCmpts;
         thisCmpts.append(word::null);
         thisCmpts.append(CompType::typeName);
         thisCmpts.append
         (
             basicThermo::splitThermoName(ThermoType::typeName(), 5)
         );
  
         wordList validNames;
  
         List<wordList> validCmpts;
         validCmpts.append
         (
             // Header
             wordList
             ({
                 typeName_(),
                 "reactionThermo",
                 "transport",
                 "thermo",
                 "equationOfState",
                 "specie",
                 "energy"
             })
         );
  
         for
         (
             const word& validName
           : dictionaryConstructorTablePtr_->sortedToc()
         )
         {
             validCmpts.append
             (
                 basicThermo::splitThermoName(validName, nCmpt)
             );
             const wordList& cmpts = validCmpts.last();
  
             bool isValid = true;
             for (label i = 1; i < cmpts.size() && isValid; ++i)
             {
                 isValid = isValid && cmpts[i] == thisCmpts[i];
             }
  
             if (isValid)
             {
                 validNames.append(cmpts[0]);
             }
         }
  
  
         FatalErrorInLookup
         (
             typeName_(),
             methodName,
             *dictionaryConstructorTablePtr_
         )
             << "All " << validCmpts[0][0] << '/' << validCmpts[0][1]
             << "/thermoPhysics combinations:" << nl << nl;
  
         // Table of available packages (as constituent parts)
         printTable(validCmpts, FatalErrorInFunction)
             << exit(FatalError);
     }
  
     return autoPtr<particleChemistryReductionMethod<CompType, ThermoType>>
     (
         cstrIter()(dict, chemistry)
     );







//    IOdictionary thermoDict
//    (
//        IOobject
//        (
//            "thermophysicalProperties",
//            dict.db().time().constant(),
//            dict.db(),
//            IOobject::MUST_READ_IF_MODIFIED,
//            IOobject::NO_WRITE,
//            false
//        )
//    );
//
//    word thermoTypeName;
//
//    if (thermoDict.isDict("thermoType"))
//    {
//        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
//        thermoTypeName =
//            word(thermoTypeDict.lookup("transport")) + '<'
//          + word(thermoTypeDict.lookup("thermo")) + '<'
//          + word(thermoTypeDict.lookup("equationOfState")) + '<'
//          + word(thermoTypeDict.lookup("specie")) + ">>,"
//          + word(thermoTypeDict.lookup("energy")) + ">";
//    }
//    else
//    {
//        FatalIOErrorInFunction(thermoDict)
//            << "thermoType is in the old format and must be upgraded"
//            << exit(FatalIOError);
//    }
//
//    dictionary MRdict(dict.subDict("reduction"));
//
//    word chemistryReductionMethodTypeName =
//        word(MRdict.lookup("method")) + '<'
//      + word(dict.subDict("chemistryType").lookup("chemistryThermo")) + ','
//      + thermoTypeName + '>';
//
//    typename dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(chemistryReductionMethodTypeName);
//
//    if (cstrIter == dictionaryConstructorTablePtr_->end())
//    {
//        FatalErrorInFunction
//            << "Unknown chemistryReductionMethodType type "
//            << chemistryReductionMethodTypeName
//            << endl << endl
//            << "Valid chemistryReductionMethodType types are :" << endl
//            << dictionaryConstructorTablePtr_->toc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<particleChemistryReductionMethod<CompType, ThermoType>>
//    (
//        cstrIter()(dict, chemistry)
//    );
}


// ************************************************************************* //
