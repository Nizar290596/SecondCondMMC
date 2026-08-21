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
    Foam::OpticallyThinRadiation

Description
    Model for emitting gas and soot species. No absorption.
    
    Gas species parameters set according to the baseline optically thin 
    radiation model recommended
    by the TNF workshop:
    http://www.sandia.gov/TNF/radiation.html
    
    Soot parameters set according to the baseline optically thin radiation 
    model recommended by the 
    International Sooting Flames Workshop
    http://www.adelaide.edu.au/cet/isfworkshop/docs/isf_baseline_radiation_model.pdf

SourceFiles
    OpticallyThinRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef OpticallyThinRadiation_H
#define OpticallyThinRadiation_H

#include "RadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class OpticallyThinRadiation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class OpticallyThinRadiation
:
    public CloudRadiationModel<CloudType>
{
    // Private data
        // Reference/Ambient Temperature
        scalar TRef_;
 
        //- Radiative species

public:

    //- Runtime type information
    TypeName("OpticallyThinRadiation");


    // Constructors

        //- Construct from dictionary
        OpticallyThinRadiation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        OpticallyThinRadiation(const OpticallyThinRadiation<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudRadiationModel<CloudType> > clone() const
        {
            return autoPtr<CloudRadiationModel<CloudType> >
            (
                new OpticallyThinRadiation<CloudType>(*this)
            );
        }


    //- Destructor
        virtual ~OpticallyThinRadiation();

    //- Member Functions

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

#ifdef NoRepository
#   include "OpticallyThinRadiation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
