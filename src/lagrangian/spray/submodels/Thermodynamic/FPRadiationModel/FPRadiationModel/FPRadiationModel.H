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
    Foam::FPRadiationModel

Description
    Templated radiation model class.

SourceFiles
    FPRadiationModel.C
    FPRadiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef FPRadiationModel_H
#define FPRadiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "vector.H"
#include "TimeFunction1.H"
#include "PtrList.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FPRadiationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FPRadiationModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("radiationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FPRadiationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        FPRadiationModel(CloudType& owner);

        //- Construct from dictionary
        FPRadiationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        FPRadiationModel(const FPRadiationModel<CloudType>& radm);

        //- Construct and return a clone
        virtual autoPtr<FPRadiationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~FPRadiationModel();


    // Selectors

        //- Selector with lookup from dictionary
        static autoPtr<FPRadiationModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );

    // Member Functions

        //- Specific Radiation Heat Flux 
        virtual scalar qR
        (
            const scalar TD,   // Droplet Temperature
            const scalar mFlux // Evaporative mass flux
        ) const = 0;
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFPRadiationModel(CloudType)                                       \
                                                                              \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                 \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        Foam::FPRadiationModel<thermoCloudType>,                              \
        0                                                                     \
    );                                                                        \
    namespace Foam                                                            \
    {                                                                         \
        defineTemplateRunTimeSelectionTable                                   \
        (                                                                     \
            FPRadiationModel<thermoCloudType>,                                \
            dictionary                                                        \
        );                                                                    \
    }


#define makeFPRadiationModelType(SS, CloudType)                               \
                                                                              \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);        \
                                                                              \
    Foam::FPRadiationModel<thermoCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType>>            \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FPRadiationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
