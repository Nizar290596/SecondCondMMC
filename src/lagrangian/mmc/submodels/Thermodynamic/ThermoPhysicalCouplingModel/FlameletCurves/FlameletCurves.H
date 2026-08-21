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
    Foam::FlameletCurves

Description
    Templated parcel single phase, multi-component class

SourceFiles
    FlameletCurves.C

\*---------------------------------------------------------------------------*/

#ifndef FlameletCurves_H
#define FlameletCurves_H

#include "ThermoPhysicalCouplingModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class FlameletCurves Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FlameletCurves
:
    public ThermoPhysicalCouplingModel<CloudType>
{
    // Private data

        //- Name of the primary species
        word primarySpecies_;
        
        //- Switch to select the independentYEqv curve
        Switch independentYEqvCurve_;
        
        //- List of TEqv target values
        List<label> TEqvTargetCurve_;
        
        //- List of YEqv target values
        List<label> YEqvTargetCurve_;
        
    // Private member functions 
        
        //- Find equivalent enthalpy and species curves
        void findEquivalentTempSpeciesCurves();

public:

    //- Runtime type information
    TypeName("FlameletCurves");

    //- Type of particle the cloud was instantiated for
    typedef typename CloudType::particleType particleType;
    
    // Constructors

        //- Construct from dictionary
        FlameletCurves
        (
            const dictionary& dict, 
            CloudType& owner, 
            const mmcVarSet& Xi
        );

    //- Destructor
    virtual ~FlameletCurves();


    //- Sources

        //- Equivalent temperature and species target values
        virtual void EqvETargetValues
        (
            const mmcVarSet& Xi,
            PtrList<volScalarField>& YEqvETarget,
            volScalarField& TEqvETarget
        );
  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FlameletCurves.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
