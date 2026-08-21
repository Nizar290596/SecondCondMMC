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
    Foam::ThermoPhysicalCouplingModel

Description
    Templated density coupling class (this allows to have different methods to
    compute target values)

SourceFiles
    ThermoPhysicalCouplingModel.C
    ThermoPhysicalCouplingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoPhysicalCouplingModel_H
#define ThermoPhysicalCouplingModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"

#include "thermoPhysicsTypes.H"
#include "interpolationLookUpTable.H"

#include "mmcVarSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ThermoPhysicalCouplingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermoPhysicalCouplingModel
:
    public SubModelBase<CloudType>
{
    // Private data

        //- Conditioning Variable in density coupling used to evaluate
        const mmcVarSet& XiC_;

        //- Coupling Names
        wordList XiCNames_;

        //- Indicator is zero or one to ensure source term is zero
        volScalarField Indicator_;

        //- Conditioning variable
        word cVarName_;

        //- List of species to solve 
        boolList solveEqvSpecie_;
        
        //- Relaxation time
        Switch tauRelaxBlending_;
        
        dimensionedScalar tauRelaxStart_;
        
        dimensionedScalar tauRelaxTarget_;
        
        scalar tauRelaxDelta_;
        
        word tauUnits_;
        
        //- Initialise thermophysical coupling data structures
        void setThermoPhysicalCoupling();

    public:

        //- Runtime type information
        TypeName("thermoPhysicalCouplingModel");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ThermoPhysicalCouplingModel,
            dictionary,
            (
                const dictionary& dict,
                CloudType& owner,
                const mmcVarSet& Xi
            ),
            (dict, owner, Xi)
        );


        // Constructors & Desctructors 

            //- Construct from owner
            ThermoPhysicalCouplingModel
            (
                CloudType& owner,
                const mmcVarSet& Xi
            );

            //- Construct from dictionary
            ThermoPhysicalCouplingModel
            (
                const dictionary& dict,
                CloudType& owner,
                const word& type,
                const mmcVarSet& Xi
            );

            //- Construct copy
            ThermoPhysicalCouplingModel(const ThermoPhysicalCouplingModel<CloudType>& cm);

            //- Construct and return a clone
            virtual autoPtr<ThermoPhysicalCouplingModel<CloudType> > clone() const
            {
                return autoPtr<ThermoPhysicalCouplingModel<CloudType> >
                       (
                            new ThermoPhysicalCouplingModel<CloudType>(*this)
                       );
            }


            //- Destructor
            virtual ~ThermoPhysicalCouplingModel();


            //- Selector
            static autoPtr<ThermoPhysicalCouplingModel<CloudType> > New
            (
                const dictionary& dict,
                CloudType& owner,
                const mmcVarSet& Xi
            );


        // Member Functions

            // Access/Edit Indicator variable
            inline const volScalarField& Indicator() const;
            inline volScalarField& Indicator();
            
        // Coupling Variables

            inline const mmcVarSet& XiC() const;

            inline const scalarField XiC0(label patch, label patchFace);
            inline const scalarField XiC0(label celli);

            //- Return the list of coupling variables name
            inline const wordList& XiCNames() const;

            //- Return the conditioning variable name
            inline const word& cVarName() const;
            
            //- Return const access to the species solved
            inline const Foam::boolList&  solveEqvSpecie() const;
            
            inline const Switch& tauRelaxBlending() const;
            
            inline const dimensionedScalar& tauRelaxStart() const;
            
            inline const dimensionedScalar& tauRelaxTarget() const;
            
            inline const scalar& tauRelaxDelta() const;
            
            inline const word& tauUnits() const;

            //- Is density coupling enabled - always returns true 
            inline bool enabled() const;

        // Sources

            //- Equivalent temperature and species target values
            //  Function is defined by the selected model
            virtual void EqvETargetValues
            (
                const mmcVarSet& Xi,
                PtrList<volScalarField>& YEqvETarget,
                volScalarField& TEqvETarget
            ) {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoPhysicalCouplingModel(CloudType)                          \
                                                                            \
    typedef CloudType::ThermoPopeCloudType ThermoPopeCloudType;             \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        ThermoPhysicalCouplingModel<ThermoPopeCloudType>,                   \
        0                                                                   \
    );                                                                      \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        ThermoPhysicalCouplingModel<ThermoPopeCloudType>,                   \
        dictionary                                                          \
    );


#define makeThermoPhysicalCouplingModelType(SS, CloudType)                  \
                                                                            \
    typedef CloudType::ThermoPopeCloudType ThermoPopeCloudType;             \
    defineNamedTemplateTypeNameAndDebug(SS<ThermoPopeCloudType>, 0);        \
                                                                            \
    ThermoPhysicalCouplingModel<ThermoPopeCloudType>::                      \
    adddictionaryConstructorToTable<SS<ThermoPopeCloudType> >               \
    add##SS##CloudType##ThermoPopeCloudType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoPhysicalCouplingModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoPhysicalCouplingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
