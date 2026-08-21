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
    Foam::CloudRadiationModel

Description
    Radiation model, coupled with soot and gas phase 

SourceFiles
    CloudRadiationModel.C
    CloudRadiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef RadiationModel_H
#define RadiationModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"

//#include "multiCmptProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CloudRadiationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudRadiationModel
:
    public SubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("radiationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CloudRadiationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CloudRadiationModel(CloudType& owner);

        //- Construct from dictionary
        CloudRadiationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        CloudRadiationModel(const CloudRadiationModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudRadiationModel<CloudType> > clone() const
        {
            return autoPtr<CloudRadiationModel<CloudType> >
            (
                new CloudRadiationModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CloudRadiationModel();


    //- Selector
    static autoPtr<CloudRadiationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Evaluation

            //- 
            virtual scalar Qrad
            (
                const scalarField& Y,
                const scalar fVsoot,
                const scalar p,
                const scalar T
            ) ;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCloudRadiationModel(CloudType)                                    \
                                                                              \
    typedef CloudType::SootPopeCloudType SootPopeCloudType;                   \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CloudRadiationModel<SootPopeCloudType>,                               \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        CloudRadiationModel<SootPopeCloudType>,                               \
        dictionary                                                            \
    );


#define makeCloudRadiationModelType(SS, CloudType)                            \
                                                                              \
    typedef CloudType::SootPopeCloudType SootPopeCloudType;                   \
    defineNamedTemplateTypeNameAndDebug(SS<SootPopeCloudType>, 0);            \
                                                                              \
    CloudRadiationModel<SootPopeCloudType>::                                  \
        adddictionaryConstructorToTable<SS<SootPopeCloudType> >               \
            add##SS##CloudType##SootPopeCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RadiationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
