/*---------------------------------------------------------------------------*\
    888b     d888                                  d8888                   
    8888b   d8888                                 d88888                   
    88888b.d88888                                d88P888                   
    888Y88888P888  .d88b.  888  888  .d88b.     d88P 888 888  888  .d88b.  
    888 Y888P 888 d88""88b 888  888 d8P  Y8b   d88P  888 888  888 d8P  Y8b 
    888  Y8P  888 888  888 Y88  88P 88888888  d88P   888 Y88  88P 88888888 
    888   "   888 Y88..88P  Y8bd8P  Y8b.     d8888888888  Y8bd8P  Y8b.     
    888       888  "Y88P"    Y88P    "Y8888 d88P     888   Y88P    "Y8888   
-------------------------------------------------------------------------------                                                                                                                                                    
License
    This file is part of movingAverage.

    movingAverage is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    movingAverage is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with movingAverage.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::aveSphere

Description
    Provide a sphere to check if point is in sphere

SourceFiles
    aveSphere.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020>

\*---------------------------------------------------------------------------*/

#ifndef aveSphere_H 
#define aveSphere_H



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "aveSpace/aveSpace/aveSpace.H"
#include "volFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class aveSphere Declaration
\*---------------------------------------------------------------------------*/

class aveSphere
:
    public aveSpace
{
    private:
    //- Private Data
        
        //- Dimensions of the sphere
        scalar radius_;

public:

    // Constructor
       
       //- Construct with three different entries
       aveSphere(const scalar radius) : radius_(radius) {};

        //- Construct from dictionary 
        //  Length must be read. If height and width are not specified, 
        //  construct equal length box
        aveSphere(const dictionary& dict);

    // Accessor functions for member variables as const reference
        
        //- Return dimensions
        inline scalar radius() const {return radius_;}

    // Helper functions

        //- Return true if the point p is in the box with the center point
        //- pCenter
        virtual bool containsPoint(const point& pCenter, const point& p) const override;

    // I-O

        //- Print box dimensions to screen
        virtual void printDim() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
