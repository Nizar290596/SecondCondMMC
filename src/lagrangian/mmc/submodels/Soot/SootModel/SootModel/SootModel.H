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
    Foam::SootModel

Description
    Templated reacting parcel composition model class
    Consists of carrier species (via thermo package), and additional liquids
    and solids

SourceFiles
    SootModel.C
    SootModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SootModel_H
#define SootModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SootModel Declaration
\*---------------------------------------------------------------------------*/
template<class CloudType>
class SootModel
:
    public SubModelBase<CloudType>
{

 public:

    //- Runtime type information
    TypeName("sootModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SootModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );

    //- Constructors

        //- Constructor null from owner
        SootModel(CloudType& owner);

        //- Constructor fom dictionary
        SootModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Constructor copy
        SootModel(const SootModel<CloudType>& cm);

        //- Constructor and return a clone
        virtual autoPtr<SootModel<CloudType> > clone() const
        {
            return autoPtr<SootModel<CloudType> >
            (
                new SootModel<CloudType>(*this)
            );
        }

    //- Destructor
        virtual ~SootModel();

    //- Selector
        static autoPtr<SootModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );

    // Member Functions
    
        // Evaluation

            // calculate the soot variables
            virtual void calSoot
            (
                scalarField& Y,
                const scalar& T,
                scalar& ySoot,
                scalar& nSoot,
                scalar& sootVf,
                scalar& sootIm,
                const scalar& rho
            );

            virtual void calSoot
            (
                scalarField& Y,
                const scalar& T,
                scalar& ySoot,
                scalar& nSoot,
                scalar& sootVf,
                scalar& sootIm,
                const scalar& rho,
                scalarList& sootContribs
            );

            //- Return if turb. mixing is switched on
            virtual Switch turbMixing() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSootModel(CloudType)                                             \
                                                                             \
    typedef CloudType::SootPopeCloudType SootPopeCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                      \
    (                                                                        \
        SootModel<SootPopeCloudType>,                                        \
        0                                                                    \
    );                                                                       \
    defineTemplateRunTimeSelectionTable                                      \
    (                                                                        \
        SootModel<SootPopeCloudType>,                                        \
        dictionary                                                           \
    );

#define makeSootModelType(SS, CloudType)                                     \
    typedef CloudType::SootPopeCloudType SootPopeCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(SS<SootPopeCloudType>, 0);           \
                                                                             \
    SootModel<SootPopeCloudType>::                                           \
    adddictionaryConstructorToTable<SS<SootPopeCloudType> >                  \
    add##SS##CloudType##SootPopeCloudType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SootModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
