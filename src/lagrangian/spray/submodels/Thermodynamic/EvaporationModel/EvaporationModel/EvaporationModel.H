/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::EvaporationModel

Description
    Templated evaporation model class

SourceFiles
    EvaporationModel.C
    EvaporationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef EvaporationModel_H
#define EvaporationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "vector.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EvaporationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EvaporationModel
:
    public CloudSubModelBase<CloudType>
{
    // Private data

public:

    //- Runtime type information
    TypeName("evaporationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        EvaporationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        EvaporationModel(CloudType& owner);

        //- Construct from dictionary
        EvaporationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        EvaporationModel(const EvaporationModel<CloudType>& evm);

        //- Construct and return a clone
        virtual autoPtr<EvaporationModel<CloudType> > clone() const
        {
            return autoPtr<EvaporationModel<CloudType> >
            (
                new EvaporationModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~EvaporationModel();


    //- Selector
    static autoPtr<EvaporationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Access

        // Evaluation
         
            //- Surface mixture fraction
            virtual scalarField YFuelVap 
            (
               const scalarField YFuel,
               const scalar Ts,
               const scalar Pa,
               const scalar Lv,
               const scalar molWtSstate
          ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelEvaporationModel(CloudType)                                 \
                                                                              \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                 \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        Foam::EvaporationModel<thermoCloudType>,                              \
        0                                                                     \
    );                                                                        \
    namespace Foam                                                            \
    {                                                                         \
        defineTemplateRunTimeSelectionTable                                   \
        (                                                                     \
            EvaporationModel<thermoCloudType>,                                \
            dictionary                                                        \
        );                                                                    \
    }                                                                         \

#define makeParcelEvaporationModelType(SS, CloudType)                         \
                                                                              \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);        \
                                                                              \
    Foam::EvaporationModel<thermoCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType> >           \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EvaporationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
