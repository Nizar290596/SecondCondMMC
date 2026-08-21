/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::ManualInjectionMMC
 
Group
    grpLagrangianIntermediateInjectionSubModels
 
Description
    Manual injection.
 
    User specifies:
      - Total mass to inject
      - Parcel positions in file \c positionsFile
      - Initial velocity field of eulerian gas, this model take the local velocity 
        of gas phase as injected droplet's initial velocity

    Properties:
      - Parcel diameters obtained by distribution model
      - All parcels introduced at SOI
 
SourceFiles
    ManualInjectionMMC.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef ManualInjectionMMC_H
#define ManualInjectionMMC_H
 
#include "InjectionModel.H"
#include "distributionModel.H"
#include "Switch.H"
#include "globalIOFields.H"
#include "volFields.H"

 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                       Class ManualInjection Declaration
\*---------------------------------------------------------------------------*/
 
template<class CloudType>
class ManualInjectionMMC
:
    public InjectionModel<CloudType>
{
    // Private data
        // Foam::autoPtr<Foam::interpolation<Foam::vector>> Uf_interp;

        // Foam::dictionary interpolationSchemes;
        //- Name of file containing the velocity field 
        const word UFieldName_;  

        // Reference to the gas phase velocity field
        const volVectorField& UField_;
 
        //- Name of file containing positions data
        const word positionsFile_;
 
        //- Field of parcel positions
        vectorGlobalIOField positions_;
 
        //- Field of parcel diameters
        scalarList diameters_;
 
        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;
 
        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;
 
        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;
 
        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;
 
        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;
 
 
public:
 
    //- Runtime type information
    TypeName("manualInjectionMMC");
 
 
    // Constructors
 
        //- Construct from dictionary
        ManualInjectionMMC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );
 
        //- Construct copy
        ManualInjectionMMC(const ManualInjectionMMC<CloudType>& im);
 
        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ManualInjectionMMC<CloudType>(*this)
            );
        }
 
 
    //- Destructor
    virtual ~ManualInjectionMMC();
 
 
    // Member Functions
 
        //- Set injector locations when mesh is updated
        virtual void updateMesh();
 
        //- Return the end-of-injection time
        scalar timeEnd() const;
 
        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);
 
        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);
 
 
        // Injection geometry
 
            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );
 
            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel

            );
 
            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;
 
            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#ifdef NoRepository
    #include "ManualInjectionMMC.C"
#endif
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* // 