/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::CompositionModel

Description
    Templated reacting parcel composition model class
    Consists of carrier species (via thermo package), and additional liquids
    and solids

    slgThermo comes from the ThermoPopeCloud class

SourceFiles
    CompositionModel.C
    CompositionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CompositionModel_H
#define CompositionModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"
#include "SLGThermo.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CompositionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CompositionModel
:
    public SubModelBase<CloudType>
{
    // Private data

        //- ThermoType  instantiated for
        typedef gasHThermoPhysics ThermoType;

        //- Reference to the thermo database
        const typename CloudType::ReactionThermoType& thermo_;

        //- SLG thermo class with access to carrier phase
        const SLGThermo& SLGThermo_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;
    
        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;


    public:

        //- Runtime type information
        TypeName("compositionModel");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            CompositionModel,
            dictionary,
            (
                const dictionary& dict,
                CloudType& owner
            ),
            (dict, owner)
        );


        // Constructors

            //- Construct null from owner
            CompositionModel(CloudType& owner);

            //- Construct from dictionary
            CompositionModel
            (
                const dictionary& dict,
                CloudType& owner,
                const word& type
            );

            //- Construct copy
            CompositionModel(const CompositionModel<CloudType>& cm);

            //- Construct and return a clone
            virtual autoPtr<CompositionModel<CloudType> > clone() const
            {
                return autoPtr<CompositionModel<CloudType> >
                       (
                            new CompositionModel<CloudType>(*this)
                       );
            }


        //- Destructor
        virtual ~CompositionModel();


        //- Selector
        static autoPtr<CompositionModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );


        // Member Functions

        // Access
            
            //- Return the thermo database
            const typename CloudType::ReactionThermoType& thermo() const;

            //- Return the SLG thermo database
            const SLGThermo& slgThermo() const;

            //- Return the carrier components (wrapper function)
            const basicMultiComponentMixture& carrier() const;


        // Mixture properties

            //- Return the list of component names for phaseI
            virtual const wordList& componentNames() const;
  
            //- Return mixture average
            virtual const ThermoType& particleMixture(const scalarField& Yp) const;

            //- Return the list of mixture mass fractions
            virtual const scalarField YMixture0(label patch,label patchFace);
            virtual const scalarField YMixture0(label celli);
            
            //- Return mole fractions given mass fractions
            scalarField X(const scalarField& Y) const;
            
            //- Return the mixture molecular weight
            scalar mixtureMW(const scalarField& Y) const;
            
            //- Return if mole fractions should be writen 
            virtual const Switch printMoleFractionsEnabled() const;


        // Specie properties
            
            //- Return molecular weights
            scalar molWt(const label i) const;

      
        // Miscellaneous
            //- Correct the mass fractions to sum to 1
            void correctMassFractions();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCompositionModel(CloudType)                                       \
                                                                              \
  typedef CloudType::ThermoPopeCloudType ThermoPopeCloudType;                 \
  defineNamedTemplateTypeNameAndDebug                                         \
  (                                                                           \
    CompositionModel<ThermoPopeCloudType>,                                    \
    0                                                                         \
  );                                                                          \
  defineTemplateRunTimeSelectionTable                                         \
  (                                                                           \
    CompositionModel<ThermoPopeCloudType>,                                    \
    dictionary                                                                \
  );


#define makeCompositionModelType(SS, CloudType)                               \
                                                                              \
  typedef CloudType::ThermoPopeCloudType ThermoPopeCloudType;                 \
  defineNamedTemplateTypeNameAndDebug(SS<ThermoPopeCloudType>, 0);            \
                                                                              \
  CompositionModel<ThermoPopeCloudType>::                                     \
  adddictionaryConstructorToTable<SS<ThermoPopeCloudType> >                   \
  add##SS##CloudType##ThermoPopeCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompositionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
