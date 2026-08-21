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
    Foam::EnvelopeModel

Description
    Templated envelope model class

SourceFiles
    EnvelopeModel.C
    EnvelopeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef EnvelopeModel_H
#define EnvelopeModel_H

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
                       Class EnvelopeModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EnvelopeModel
:
    public CloudSubModelBase<CloudType>
{
    // Private data


public:

    //- Runtime type information
    TypeName("envelopeModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        EnvelopeModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );
protected:

public:
    // Constructors

        //- Construct null from owner
        EnvelopeModel(CloudType& owner);

        //- Construct from dictionary
        EnvelopeModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        EnvelopeModel(const EnvelopeModel<CloudType>& efm);

        //- Construct and return a clone
        virtual autoPtr<EnvelopeModel<CloudType> > clone() const
        {
            return autoPtr<EnvelopeModel<CloudType> >
            (
                new EnvelopeModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~EnvelopeModel();


    //- Selector
    static autoPtr<EnvelopeModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );

//**//    scalarField f1_;
//**//    scalarField f0_;
    scalar f1_; //**//
    scalar f0_; //**//

    scalar TIgnition_;

    // Member Functions

        // Access


	// Edit

//**//        inline scalarField& f1();
//**//        inline scalarField& f0();
        inline scalar& f1(); //**//
        inline scalar& f0(); //**//

        // Evaluation

            //- initialize parcel
            virtual void initialize
            (
               const scalarField Y0,
               const scalar T0,
	       const label sCell
            );

            virtual void setSwitches();

            virtual bool checkState
            (
               const scalar T0
	        );

            virtual void setBoundaries
            (
               const scalar f0,
               const scalar f1
            );

            //- mixture fraction: no envelope
//**//            virtual scalarField fNoEnvelope
            virtual scalar fNoEnvelope //**//
            (
               const scalarField Y0,
               const scalarField Y1,
               const scalarField YX
            );
        
            //- compostion: no envelope
            virtual scalarField YNoEnvelope
            (
               const scalarField Y0,
               const scalarField Y1,
               const scalar fX
            ) const;

            //- total enthalpy: no envelope
            virtual scalar haNoEnvelope
            (
               const scalarField YX,
               const scalar p,
               const scalar TD
            ) const;

            //- temperature: no envelope
            virtual scalar TNoEnvelope
            (
               const scalar TD
            ) const;

//**//            virtual scalarField fSurface
            virtual scalar fSurface //**//
            (
                const scalarField Y0,
                const scalarField Y1,
                const scalar f0,
                const scalarField YSurf
            );

            virtual scalarField YSurface
            (
                const scalarField Y0Liq,
                const scalarField Y1,
                const scalar fSurf
            );

            virtual scalarField YFilm
            (
                const scalarField Y0Liq,
                const scalarField Y1,
                const scalar fFilm
            ) const;

	    virtual scalar haSurface
            (
                const scalarField YSurf,
                const scalar p,
                const scalar TD,
                const scalar ha0,
                const scalar ha1,
                const scalar haSurf
	    ) const;

	    virtual scalar TSurface
            (
                const scalarField YSurf,
                const scalar p,
                const scalar Td,
                const scalar TGas,
                const scalar fSurf,
                const scalar haSurf
            ) const;

        //- turns composition of n_fuel components into n_thermo components 
	    virtual scalarField Y1
	    (
	        const scalarField YFuelVap	
	    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeEnvelopeModel(CloudType)                                        \
                                                                            \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;               \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        Foam::EnvelopeModel<thermoCloudType>,                               \
        0                                                                   \
    );                                                                      \
    namespace Foam                                                          \
    {                                                                       \
        defineTemplateRunTimeSelectionTable                                 \
        (                                                                   \
            EnvelopeModel<thermoCloudType>,                                 \
            dictionary                                                      \
        );                                                                  \
    }                                                                       \

#define makeEnvelopeModelType(SS, CloudType)                                \
                                                                            \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;               \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);      \
                                                                            \
    Foam::EnvelopeModel<thermoCloudType>::                                  \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType> >         \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EnvelopeModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
