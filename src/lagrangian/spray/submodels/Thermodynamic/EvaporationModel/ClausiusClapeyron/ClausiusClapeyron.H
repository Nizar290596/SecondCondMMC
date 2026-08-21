/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::ThermoSurfaceFilm

Description
    Thermo parcel surface film model.

    Responsible for:
    - injecting parcelss from the film model into the cloud, e.g. for dripping
    - parcel interaction with the film, e.g absorb, bounce, splash

    Splash model references:

        Bai and Gosman, `Mathematical modelling of wall films formed by
        impinging sprays', SAE 960626, 1996

        Bai et al, `Modelling off gasoline spray impingement', Atom. Sprays,
        vol 12, pp 1-27, 2002


SourceFiles
    ClausiusClapeyron.C
    ClausiusClapeyron.H

\*---------------------------------------------------------------------------*/

#ifndef ClausiusClapeyron_H
#define ClausiusClapeyron_H

#include "EvaporationModel.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThermoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClausiusClapeyron
:
   public EvaporationModel<CloudType>
{

private: 

    // privated data

        //- Molecular weight of species
        const scalar molWtDstate_;
    
        //-Boilinig temperature of fuel parcel
        const scalar Tb_;

        //- Universal gas constant
        const scalar Rg_;   
       
        //- Reference pressure for T boil
	    const scalar pRef_;
public:

        //- Runtime type information
        TypeName("ClausiusClapeyron");


        // Constructors

        //- Construct from components
        ClausiusClapeyron(const dictionary& dict, CloudType& owner);

        //- Construct copy
        ClausiusClapeyron(const ClausiusClapeyron<CloudType>& ccp);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<EvaporationModel<CloudType> > clone() const
        {
            return autoPtr<EvaporationModel<CloudType> >
            (
                new ClausiusClapeyron<CloudType>(*this)
            );
        }

        
        //- Destructor
        virtual ~ClausiusClapeyron();


        // Member Functions
    
        // Access

 
        //-Surface mixture fraction
          
        virtual scalarField YFuelVap
            (
                const scalarField YFuel,//- composition of the liquid fuel
                const scalar Ts,        //- Surface temperature
                const scalar Pa,        //- Refernce pressure for TBoil
                const scalar Lv,        //- latent heat
                const scalar molWtSstate//- molecular weight of surface state 
            ) const;
            
     
        // I-O
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#ifdef NoRepository
#   include "ClausiusClapeyron.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
