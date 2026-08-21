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

\*---------------------------------------------------------------------------*/

#include "ClausiusClapeyron.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Pstream.H"
#include "interpolationLookUpTable.H"
#include "particle.H"
#include "cloud.H"    


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ClausiusClapeyron<CloudType>::ClausiusClapeyron
(
    const dictionary& dict,
    CloudType& owner
)
:
   
EvaporationModel<CloudType>(dict,owner,typeName),
molWtDstate_(this->owner().constProps().molWtSstate0()),
Tb_(readScalar(this->coeffDict().lookup("Tb"))),
Rg_(readScalar(this->coeffDict().lookup("Rgas"))),
pRef_(readScalar(this->coeffDict().lookup("pRef")))
{}


template <class CloudType>
Foam::ClausiusClapeyron<CloudType>::ClausiusClapeyron
(
    const ClausiusClapeyron<CloudType>& cclap
)
:

EvaporationModel<CloudType>(cclap),
molWtDstate_(readScalar(this->coeffDict().lookup("molWtDstate"))),
Tb_(readScalar(this->coeffDict().lookup("Tb"))),  
Rg_(readScalar(this->coeffDict().lookup("Rgas"))),
pRef_(readScalar(this->coeffDict().lookup("pRef")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClausiusClapeyron<CloudType>::~ClausiusClapeyron()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Calculate mixture fraction at surface
template<class CloudType>
Foam::scalarField Foam::ClausiusClapeyron<CloudType>::YFuelVap
(
    const scalarField YFuel,
    const scalar Ts,
    const scalar Pa,
    const scalar Lv,
    const scalar molWtSstate
) const
{   
        scalarField YFuelSurf;
        YFuelSurf.setSize(YFuel.size(),0.0);

        forAll(YFuel,nfs)
        {
            scalar constantB = Lv * molWtDstate_ / Rg_;
            scalar constantA = pRef_ * exp( constantB / Tb_ );
            scalar P_FS =  constantA * exp( - constantB / Ts );
            YFuelSurf[nfs] = (1 / ( 1 + molWtSstate / molWtDstate_ * ( Pa / P_FS - 1 ) ));

            YFuelSurf[nfs] = min(YFuelSurf[nfs],0.99); // Y >= 1.0 can cause problems
            
            YFuelSurf[nfs] *= YFuel[nfs];              // Scales YFuelSurf in case for multicomponent fuels
                                                       // (needs to be changes for differential evaporation),
                                                       // has no effect on single component fuels
        }

        return YFuelSurf;
}
// ************************************************************************* //
