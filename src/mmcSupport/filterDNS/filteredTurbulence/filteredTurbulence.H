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
    Foam::filteredTurbulence

Description
    Generates all filtered turbulent properties.

    Automatic handling of the flag filterDNSFields. 

SourceFiles
    filteredTurbulence.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef filteredTurbulence_H
#define filteredTurbulence_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "filterDNSField.H"
#include "psiReactionThermo.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                        Class filteredTurbulence Declaration
\*---------------------------------------------------------------------------*/
class filteredTurbulence
{
    private:
     
    // Private Member Variables
        
        //- Filter model
        filterDNSField filter_;

        //- Thermo model required for turbulence model
        autoPtr<psiReactionThermo> thermo_;

        //- Refernce to the filterered species mass fraction
        PtrList<volScalarField>& Y_;
        
        //- Switch if filtered fields are used for the mmcFoam cloud or not
        const bool filteredMMCCloudFields_;

        //- Laminar and turbulent Schmidt number
        const scalar Sc_;
        const scalar Sct_;

        //- autoPtr to the filtered fields
        autoPtr<volScalarField> rhoFilteredPtr_;
        autoPtr<volVectorField> UFilteredPtr_;
        autoPtr<surfaceScalarField> phiFilteredPtr_;

        //- autoPtr to the turbulent diffusivity fields
        autoPtr<volScalarField> D_;
        autoPtr<volScalarField> Dt_;
        autoPtr<volScalarField> DEff_;
        autoPtr<volScalarField> DeltaE_;

        //- Turbulence model
        autoPtr<compressible::turbulenceModel> turbulence_;



    // Private Member Functions

        
    
    public:

        // Constructors

            filteredTurbulence
            (
                const fvMesh& mesh,
                const fvMesh& filterMesh,
                const volScalarField& rhoDNS,
                const volVectorField& UDNS,
                const bool filteredMMCCloudFields,
                const scalar Sc,
                const scalar Sct
            );

            //- Disallow default bitwise copy construct
            filteredTurbulence(const filteredTurbulence&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const filteredTurbulence&) = delete;

        // Correct 

            //- Correct the turbulence and thermo model based on the DNS 
            //  field information
            void correct
            (
                const volVectorField& UDNS,
                const volScalarField& rhoDNS,
                const volScalarField& pDNS,
                const volScalarField& TDNS,
                const PtrList<volScalarField>& Y
            );

        // Const Access

            //- Access to the turbulence model
            const compressible::turbulenceModel& turbModel() {return turbulence_();}

            //- Access to the thermo model 
            const psiReactionThermo& thermo() {return thermo_();}

            //- Access to the filter function
            const filterDNSField& filter() const {return filter_;}

            //- Access to the filtered velocity field
            //  Note this is always a filtered velocity also if filteredMMCCloudFields
            //  is switched off
            const volVectorField& U() {return UFilteredPtr_();}

            //- Access to the filtered density field
            //  Note this is always a filtered density also if filteredMMCCloudFields
            //  is switched off
            const volScalarField& rho() {return rhoFilteredPtr_();}

            //- Access to the turbulent diffusivity fields
            const volScalarField& D() {return D_();}
            const volScalarField& Dt() {return Dt_();}
            const volScalarField& DEff() {return DEff_();}
};

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
