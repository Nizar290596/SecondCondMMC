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
    Foam::nameVariableTable
 
Description
    Class provides a function to access variables stored as references by 
    a given name. E.g. storing a variable rho with the name "density" to 
    have a descriptive access. 

    This is used for the particle statistics function, allowing users to 
    select which variables to sample.

SourceFiles
    nameVariableTable.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024

\*---------------------------------------------------------------------------*/

#ifndef nameVariableTable_H
#define nameVariableTable_H

#include "fvCFD.H"
#include <vector>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class nameVariableTable
\*---------------------------------------------------------------------------*/
// Base struct to store the eulerian field properties for each particle

class nameVariableTable
{
    private:
        
        //- Static variable to point to the list of hash tables containing 
        //  the HashTable for the current nameVariableTable.
        //  Key is the name of the nameVariableTable stored in tableName_
        //  The concept is, that the names are stored only once as they
        //  use a lot of space. 
        static HashTable<HashTable<size_t>> nameIndexMapTable_;

        //- Name of this nameVariableTable
        const word tableName_{"default"};

        //- Store the variable names as a static HashTable
        //  The concept is, that the names are stored only once as they
        //  use a lot of space. 
        //  In the map the variable name is associated to an index for the 
        //  lookup in the std::vector<scalar*> varPointers_
        HashTable<size_t>& nameIndexMap_;

        //- Storing the pointers to the variables
        std::vector<const scalar*> varPointers_;

        //- Get the nameIndexMap from nameIndexMapTable
        //  Only used during constructor call
        HashTable<size_t>& getNameIndexMap();

        //- Check if var name exists in the nameIndexMap.
        //  If not it is added and the index is returned
        size_t getNameOrAddIndex(const word& varName);

        //- Check if var name exists in the nameIndexMap.
        size_t getNameIndex(const word& varName) const;
    
    public:
    
    // Constructors

        //- Default constructor
        nameVariableTable();

        //- Construct with given name
        nameVariableTable(const word tableName);

    // Public Member Functions

        //- Get the name of this nameVariableTable
        //  Typically uses the name of the Lagrangian parcel that is tracked
        const word& name() const {return tableName_;}

        //- Add a new variable
        void addNamedVariable(const word& varName, const scalar& var);

        //- Get access to the stored variable
        const scalar& get(const word& varName) const;

        //- Get the names of all variables stored
        wordList getAllVarNames() const;

        //- Add all variables of the word list to a container
        //  The container type must provide an Type::append() function
        template<class Type>
        void storeVarsByList(Type& container, const wordList& vars) const;

        //- Add all variables stored to the storage container
        //  The container type must provide an Type::append() function
        template<class Type>
        void storeAllVars(Type& container) const;
};

}   // End of namespace Foam

template<class Type>
void Foam::nameVariableTable::storeVarsByList
(
    Type& container,
    const wordList& vars
) const
{
    container.reserve(container.size()+vars.size());
    for (const word& var : vars)
    {
        auto ind = nameVariableTable::getNameIndex(var);
        container.append(*(varPointers_[ind]));
    }
}


template<class Type>
void Foam::nameVariableTable::storeAllVars
(
    Type& container
) const
{
    container.reserve(container.size()+varPointers_.size());
    for (auto it=nameIndexMap_.begin(); it != nameIndexMap_.end();it++)
    {
        container.append(*(varPointers_[*it]));
    }
}






#endif
