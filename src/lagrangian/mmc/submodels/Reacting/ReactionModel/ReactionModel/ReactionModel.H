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
    Foam::ReactionModel

Description
    Templated reacting parcel composition model class
    Consists of carrier species (via thermo package), and additional liquids
    and solids

SourceFiles
    ReactionModel.C
    ReactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ReactionModel_H
#define ReactionModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ReactionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactionModel
:
    public SubModelBase<CloudType>
{
    // Private data

public:

    //- Runtime type information
    TypeName("reactionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ReactionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        ReactionModel(CloudType& owner);

        //- Construct from dictionary
        ReactionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        ReactionModel(const ReactionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<ReactionModel<CloudType> > clone() const
        {
            return autoPtr<ReactionModel<CloudType> >
            (
                new ReactionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ReactionModel() = default;


    //- Selector
    static autoPtr<ReactionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Ignite the flamesheet
        virtual void ignite(){}
        
        //- Calcuate chemistry
        virtual void calculate
        (
            const scalar t0, 
            const scalar dt,
            const scalar ha, 
            const scalar pc, 
            const scalar pT,
            const scalar pz,
            scalarField& Y
        ){}

        virtual void ISATCntIncrmnt(){}

        virtual void resetCPUTimeAnalysis(){}

        virtual void writeCPUTimeAnalysis(){}

        virtual label nActiveSp(){ return 0;}

        virtual void Sreact(){}

        //virtual scalar getCPUTimeReaction(){ return 0.0;}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactionModel(CloudType)                                        \
                                                                            \
    typedef CloudType::reactingPopeCloudType reactingPopeCloudType;         \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        ReactionModel<reactingPopeCloudType>,                               \
        0                                                                   \
    );                                                                      \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        ReactionModel<reactingPopeCloudType>,                               \
        dictionary                                                          \
    );


#define makeReactionModelType(SS,CloudType)                                 \
                                                                            \
    typedef CloudType::reactingPopeCloudType reactingPopeCloudType;         \
    defineNamedTemplateTypeNameAndDebug(SS<reactingPopeCloudType>, 0);      \
                                                                            \
    ReactionModel<reactingPopeCloudType>::                                  \
      adddictionaryConstructorToTable<SS<reactingPopeCloudType>>            \
            add##SS##CloudType##reactingPopeCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
