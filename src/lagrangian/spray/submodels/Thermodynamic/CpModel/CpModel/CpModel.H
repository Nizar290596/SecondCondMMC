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
    Foam::CpModel

Description
    Templated surface enthalpy calculation class.



\*---------------------------------------------------------------------------*/

#ifndef CpModel_H
#define CpModel_H

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
                       Class CpModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CpModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("CpModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CpModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CpModel(CloudType& owner);

        //- Construct from dictionary
        CpModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        CpModel(const CpModel<CloudType>& cpm);

        //- Construct and return a clone
        virtual autoPtr<CpModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~CpModel();


    // Selectors

        //- Selector with lookup from dictionary
        static autoPtr<CpModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );

    // Member Functions

        //- Specific Heat at Constant Pressure 
        virtual scalar Cp
        (
            const scalar T
        ) const = 0;
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCpModel(CloudType)                                                 \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CpModel<thermoCloudType>,                                        \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CpModel<thermoCloudType>,                                          \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCpModelType(SS, CloudType)                                         \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);         \
                                                                               \
    Foam::CpModel<thermoCloudType>::                                           \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType>>             \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CpModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
