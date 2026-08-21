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
    Foam::particleChemistryModel
 
Description
    Extends the StandardChemistryModel of OpenFOAM with a particleCalculate
    function that solves the chemistry on a particle level. 

    Note:
        In old mmcFoam this was called popeParticleChemistryModel

SourceFiles
    particleChemistryModel.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/
 
#ifndef particleChemistryModel_H
#define particleChemistryModel_H
 
#include "StandardChemistryModel.H"
#include "OFstream.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                      Class particleChemistryModel Declaration
\*---------------------------------------------------------------------------*/
 
template<class ReactionThermo, class ThermoType>
class particleChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{

public:
    //- Runtime type information
    TypeName("particleChemistryModel");

    // Constructors

        //- No copy construct
        particleChemistryModel(const particleChemistryModel&) = delete;
  
        //- No copy assignment
        void operator=(const particleChemistryModel&) = delete;
  
        //- Construct from thermo
        particleChemistryModel(ReactionThermo& thermo);
  
  
        //- Destructor
        virtual ~particleChemistryModel() {};
        
        //- Calculate the estimation of chemical integration time step
        //  Note: pi and Ti cannot be constant because of the call in the 
        //        ode solver -- see ode.H
        virtual scalar particleCalculate
        (
            const scalar t0,        // initial time
            const scalar deltaT,    // time step
            const scalar hi,        // enthalpy - currently not used
            scalar pi,              // pressure of particle i
            scalar Ti,              // temperature of particle i
            scalarField& Y          // list of species of particle i
        );
};

}   // End of namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleChemistryModel.C"
#endif

#endif



