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

\*---------------------------------------------------------------------------*/

#include "nameVariableTable.H"


// * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * * * *
Foam::HashTable<HashTable<size_t>> Foam::nameVariableTable::nameIndexMapTable_;


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * *

Foam::nameVariableTable::nameVariableTable()
:
nameIndexMap_(getNameIndexMap())
{

}


Foam::nameVariableTable::nameVariableTable(const word tableName)
:
tableName_(tableName),
nameIndexMap_(getNameIndexMap())
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * *

Foam::HashTable<size_t>& Foam::nameVariableTable::getNameIndexMap()
{
    auto it = nameIndexMapTable_.find(tableName_);
    if (it != nameIndexMapTable_.end())
    {
        return *it;
    }
    else
    {
        nameIndexMapTable_.insert(tableName_,HashTable<size_t>());
        return *(nameIndexMapTable_.find(tableName_));
    }

    // Dummy return to avoid warnings during compilation
    return *it;
}



size_t Foam::nameVariableTable::getNameOrAddIndex(const word& varName)
{
    const auto nVars = varPointers_.size();
    auto it = nameIndexMap_.find(varName);
    if (it != nameIndexMap_.end())
        return *it;
    else
    {
        auto res = nameIndexMap_.insert(varName,nVars);
        if (res)
            return nVars;
    }

    FatalError << "Could not insert varName: "<<varName 
        << " in Foam::nameVariableTable::getNameIndex()"
        << exit(FatalError);

    return 0;
}


size_t Foam::nameVariableTable::getNameIndex(const word& varName) const
{
    auto it = nameIndexMap_.find(varName);
    
    #ifdef FULLDEBUG
    if (it == nameIndexMap_.end())
    {
        List<word> toc = nameIndexMap_.toc();
        FatalError << "Variable " << varName
                   << " is not available on the particle and nameVariableTable "
                   << tableName_ << nl
                   << "Following variables can be selected:" << nl
                   << toc << exit(FatalError);
    }
    #endif
    return *it;
}


void Foam::nameVariableTable::addNamedVariable
(
    const word& varName,
    const scalar& var
)
{
    // Get the index in the vector for this name
    auto ind = nameVariableTable::getNameOrAddIndex(varName);

    // Check if index is accessible 
    if (ind >= varPointers_.size())
    {
        varPointers_.reserve(ind+1);
        varPointers_.resize(ind+1,nullptr);
    }

    varPointers_[ind] = &var;
}


const scalar& Foam::nameVariableTable::get(const word& varName) const
{
    auto ind = getNameIndex(varName);
    #ifdef FULLDEBUG
    if (ind >= varPointers_.size())
    {
        FatalError << "Mismatch between varPointers and nameVariableTable for "
                   << "variable " << varName << nl
                   << nameIndexMap_.size() << " variable names are stored "
                   << "but only " << varPointers_.size() << " variable entries"
                   << exit(FatalError);
    }
    #endif

    return *(varPointers_[ind]);
}


Foam::wordList Foam::nameVariableTable::getAllVarNames() const
{
    List<word> varNames(nameIndexMap_.size());
    label i=0;
    for (auto it=nameIndexMap_.begin(); it != nameIndexMap_.end();it++)
    {
        varNames[i++] = it.key();
    }
    
    return varNames;
}

