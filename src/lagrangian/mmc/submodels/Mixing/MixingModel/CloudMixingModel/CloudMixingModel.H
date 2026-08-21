/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::MixingModel

Description
    Templated Mixing model class

SourceFiles
    MixingModel.C
    MixingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef MixingModel_H
#define MixingModel_H

#include "SubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"
#include "mmcVarSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MixingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudMixingModel
:
    public SubModelBase<CloudType>
{
    public:
    
        //- Runtime type information
        TypeName("mixingModel");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            CloudMixingModel,
            dictionary,
            (
                const dictionary& dict,
                CloudType& owner,
                const mmcVarSet& Xi //!!!
            ),
            (dict, owner, Xi)
        );


        // Constructors

            //- Construct null from owner
            CloudMixingModel(CloudType& owner); // !!!

            //- Construct from dictionary
            CloudMixingModel
            (
                const dictionary& dict,
                CloudType& owner,
                const word& type
            );

            //- Construct copy
            CloudMixingModel(const CloudMixingModel<CloudType>& cm);

            //- Construct and return a clone
            virtual autoPtr<CloudMixingModel<CloudType> > clone() const
            {
                return autoPtr<CloudMixingModel<CloudType> >
                (
                    new CloudMixingModel<CloudType>(*this)
                );
            }


        // Destructor
        
            //- Default destructor 
            virtual ~CloudMixingModel() {};


        //- Selector
        static autoPtr<CloudMixingModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi //!!!
        );

        // Member Functions

            //- Return mixing impulse 
            //- scalar values updated along with integration of chemical source term
            virtual void Smix();
            
            //- Return scalar field with reference variables
            virtual inline const scalarField XiR0(label patch, label patchFace, particle& p);
            virtual inline const scalarField XiR0(label celli, particle& p);
            
            virtual inline const wordList& XiRNames() const;
            
            //- Number of reference variables 
            virtual inline label numXiR() const;
            
            //- return mmcVarSet
            virtual const mmcVarSet& XiR() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCloudMixingModel(CloudType)                                     \
                                                                            \
    typedef CloudType::mixingPopeCloudType mixingPopeCloudType;             \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        CloudMixingModel<mixingPopeCloudType>,                              \
        0                                                                   \
    );                                                                      \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        CloudMixingModel<mixingPopeCloudType>,                              \
        dictionary                                                          \
    );


#define makeCloudMixingModelType(SS, CloudType)                             \
                                                                            \
    typedef CloudType::mixingPopeCloudType mixingPopeCloudType;             \
    defineNamedTemplateTypeNameAndDebug(SS<mixingPopeCloudType>, 0);        \
                                                                            \
    CloudMixingModel<mixingPopeCloudType>::                                 \
        adddictionaryConstructorToTable<SS<mixingPopeCloudType> >           \
            add##SS##CloudType##mixingPopeCloudType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudMixingModel.C"
    #include "CloudMixingModelNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
