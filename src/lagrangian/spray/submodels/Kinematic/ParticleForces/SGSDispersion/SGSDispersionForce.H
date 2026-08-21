/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::SGSDispersionForce

Description
    Calculates stochastic particle dispersion component based on the model of 
    Bini and Jones.
    

    Reference:
    \verbatim
        [1] Bini and Jones (2007), Particle acceleration in turbulent flows: A 
            class of nonlinear stochastic models for intermittency.
        [2] Bini and Jones (2008), Large-eddy simulation of particle-laden turbulent flows.   
        [3] Bini and Jones (2009), Large Eddy Simulation of an evaporating 
            acetone spray.
    \endverbatim

SourceFiles
    SGSDispersionForceI.H
    SGSDispersionForce.C

\*---------------------------------------------------------------------------*/

#ifndef SGSDispersionForce_H
#define SGSDispersionForce_H

#include "ParticleForce.H"
#include "Random.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SGSDispersionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SGSDispersionForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Reference to the cloud random number generator
        Random& rndGen_;

        //- Dispersion Constant []
        const scalar C0_;

        //- Turbulent time scale exponent parameter
        const scalar alpha_;

        //- Turbulent time scale time scale type
        const word tauTName_;

        //- LES Delta
        const volScalarField& Delta_; 

        //- Pointer to the turbulence kinetic energy field
        const volScalarField* kPtr_;

        //- Flag that indicates ownership of turbulence k field
        bool ownK_;


    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;

        //- Weiner process increment
        vector dW(const scalar dt) const;

        //- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;


public:

    //- Runtime type information
    TypeName("sgsDispersion");


    // Constructors

        //- Construct from mesh
        SGSDispersionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        SGSDispersionForce(const SGSDispersionForce& bmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new SGSDispersionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SGSDispersionForce();


    // Member Functions

        // Access

            //- Return const access to the Dispersion Constant []
            inline scalar C0() const;

            //- Return const access to Turbulent time scale exponent parameter []
            inline scalar alpha() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SGSDispersionForceI.H"

#ifdef NoRepository
    #include "SGSDispersionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
