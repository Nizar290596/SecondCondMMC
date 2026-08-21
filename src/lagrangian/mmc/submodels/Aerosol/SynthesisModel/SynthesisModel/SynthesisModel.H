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
    Foam::SynthesisModel

Description
    Templated synthesis model class. 

SourceFiles
    SynthesisModel.C
    SynthesisModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SynthesisModel_H
#define SynthesisModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SynthesisModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SynthesisModel
:
    public SubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("synthesisModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SynthesisModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        SynthesisModel(CloudType& owner);

        //- Construct from dictionary
        SynthesisModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        SynthesisModel(const SynthesisModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SynthesisModel<CloudType> > clone() const
        {
            return autoPtr<SynthesisModel<CloudType> >
            (
                new SynthesisModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SynthesisModel();


    //- Selector
    static autoPtr<SynthesisModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // synthesis
        virtual void synthesize();
        
        virtual const wordList psdPropertyNames() const;
        
        virtual const wordList physicalAerosolPropertyNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSynthesisModel(CloudType)                                        \
                                                                             \
    typedef CloudType::aerosolPopeCloudType aerosolPopeCloudType;        \
    defineNamedTemplateTypeNameAndDebug                                      \
    (                                                                        \
        SynthesisModel<aerosolPopeCloudType>,                              \
        0                                                                    \
    );                                                                       \
    defineTemplateRunTimeSelectionTable                                      \
    (                                                                        \
        SynthesisModel<aerosolPopeCloudType>,                              \
        dictionary                                                           \
    );

#define makeSynthesisModelType(SS, CloudType)                                \
                                                                             \
    typedef CloudType::aerosolPopeCloudType aerosolPopeCloudType;        \
    defineNamedTemplateTypeNameAndDebug(SS<aerosolPopeCloudType>, 0);      \
                                                                             \
    SynthesisModel<aerosolPopeCloudType>::                                 \
        adddictionaryConstructorToTable<SS<aerosolPopeCloudType> >         \
            add##SS##CloudType##aerosolPopeCloudType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SynthesisModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
