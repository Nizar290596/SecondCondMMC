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
    Foam::baseParticleDataContainer
 
Description
    Base class to store particle data in a new container that can be used for 
    easier inter processor streaming with Pstream.
     
    The data container contains two main lists, one for scalars, the other for
    scalarFields. The number of scalars and scalarFields has to be set when
    the main particleDataContainer is derived from this base class. 
    
    See also mixParticleModel.H

SourceFiles
    baseParticleDataContainer.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/

#ifndef baseParticleDataContainer_H
#define baseParticleDataContainer_H

#include "fvCFD.H"
#include "Ostream.H"
#include "Istream.H"


namespace Foam
{
/*---------------------------------------------------------------------------*\
                IO Functions required for Parallel Handling
\*---------------------------------------------------------------------------*/

class baseParticleDataContainer;

// Read and write the eulerian fields data class
Istream& operator >> (Istream& is, baseParticleDataContainer& eField);

Ostream& operator << (Ostream& os, const baseParticleDataContainer& eField);

// When comparing baseParticleDataContainer always return untrue as otherwise all 
// fields have to be checked
// --> This operation is only needed to write out a list of baseParticleDataContainer
//     which checks for uniform list entries
inline bool operator==(const baseParticleDataContainer& e1, const baseParticleDataContainer& e2)
{
    return false;
}

inline bool operator!=(const baseParticleDataContainer& e1, const baseParticleDataContainer& e2)
{
    return true;
}


/*---------------------------------------------------------------------------*\
                     Class baseParticleDataContainer
\*---------------------------------------------------------------------------*/
// Base struct to store the eulerian field properties for each particle

class baseParticleDataContainer
{
    protected:
        
        //- General list of scalar values to store
        List<label> labelValues_;
        
        //- General list of scalar values to store
        List<scalar> scalarValues_;
        
        //- General list of scalarLists or scalarFields to store
        List<Field<scalar>> scalarFields_;
        
        //- Position of the particle in the x,y,z space
        vector pos_{0,0,0}; 

        //- Is it a local particle
        //  set to false if transferred
        bool local_{true};
    
    public:
    
        //- Construct with given size for lists 
        baseParticleDataContainer
        (
            const label nLabelValues,
            const label nScalarValues,
            const label nScalarLists
        )
        : 
            labelValues_(nLabelValues), 
            scalarValues_(nScalarValues), 
            scalarFields_(nScalarLists) 
        {};
        
        //- Construct from input stream
        baseParticleDataContainer(Istream& is);
        
        // Functions for IO operations
        
            //- Read baseParticleDataContainer from stream
            Istream& read(Istream& is);
            
            //- Write baseParticleDataContainer to stream
            Ostream& write(Ostream& os) const;
           

       // Modify

            //- Access stored content directly
            List<label>& labelValues() {return labelValues_;}
            List<scalar>& scalarValues() {return scalarValues_;}
            List<Field<scalar>>& scalarFields() {return scalarFields_;}
            
            //- Const access stored content directly
            const List<label>& labelValues() const {return labelValues_;}
            const List<scalar>& scalarValues() const {return scalarValues_;}
            const List<Field<scalar>>& scalarFields() const {return scalarFields_;}



            //- Add two baseParticleDataContainers
            //  Position variable is unmodified
            void append(baseParticleDataContainer& container);

        // Edit
        
            //- Non const access to the position
            inline vector& position() {return pos_;}

            inline bool& local() {return local_;}
        
        // Access

            //- Const Access to the position
            inline const vector& position() const {return pos_;}

            //- Const access if it is a local particle
            inline const bool& local() const {return local_;}
};

}   // End of namespace Foam
#endif
