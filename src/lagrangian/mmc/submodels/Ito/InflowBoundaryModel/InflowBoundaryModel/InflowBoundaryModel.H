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
    Foam::InflowBoundaryModel


Description
    Templated inflow boundary model class

SourceFiles
    InflowBoundaryModel.C
    NewInflowBoundaryModel.C

\*---------------------------------------------------------------------------*/

#ifndef InflowBoundaryModel_H
#define InflowBoundaryModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
//#include "addToRunTimeSelectionTable.H"
#include "SubModelBase.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class InflowBoundaryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InflowBoundaryModel
  :
    public SubModelBase<CloudType>
{

public:
   //- Convenience typedef for parcelType
    typedef typename CloudType::particleType particleType;

    //- Runtime type information
    TypeName("InflowBoundaryModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InflowBoundaryModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        InflowBoundaryModel(CloudType& owner);

        //- Construct from dictionary
        InflowBoundaryModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );
        //- Construct copy
        InflowBoundaryModel(const InflowBoundaryModel<CloudType>& ibc);

        //- Construct and return a clone
        virtual autoPtr<InflowBoundaryModel<CloudType> > clone() const
        {
            return autoPtr<InflowBoundaryModel<CloudType> >
            (
                new InflowBoundaryModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InflowBoundaryModel();


    //- Selector
    static autoPtr<InflowBoundaryModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

    virtual void inflow();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInflowBoundaryModel(CloudType)                                  \
                                                                            \
    typedef CloudType::itoPopeCloudType itoPopeCloudType;                   \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        InflowBoundaryModel<itoPopeCloudType>,                              \
        0                                                                   \
    );                                                                      \
                                                                            \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        InflowBoundaryModel<itoPopeCloudType>,                              \
        dictionary                                                          \
    );


#define makeInflowBoundaryModelType(SS, CloudType)                          \
                                                                            \
    typedef CloudType::itoPopeCloudType itoPopeCloudType;                   \
    defineNamedTemplateTypeNameAndDebug(SS<itoPopeCloudType>, 0);           \
                                                                            \
    InflowBoundaryModel<itoPopeCloudType>::                                 \
        adddictionaryConstructorToTable<SS<itoPopeCloudType> >              \
            add##SS##CloudType##itoPopeCloudType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InflowBoundaryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
