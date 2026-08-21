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

\*---------------------------------------------------------------------------*/

#include "newtonLinODEExtern.H"
#include "finiteRateNewtonLinODE.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::finiteRateNewtonLinODE<CloudType>::finiteRateNewtonLinODE
(
    const dictionary& dict,
    CloudType& owner
)
:
    ode
    <
        particleChemistryModel
        <
            psiReactionThermo,
            gasHThermoPhysics
        >
    >(dynamic_cast<psiReactionThermo&>(owner.thermo())),
    
    ReactionModel<CloudType>(dict,owner,typeName),

    zLower_(this->coeffDict().lookupOrDefault("zLower", 0.0)),
    
    zUpper_(this->coeffDict().lookupOrDefault("zUpper", 1.0))
{
    newtonlinodeinit_
    (
    );
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::finiteRateNewtonLinODE<CloudType>::~finiteRateNewtonLinODE()
{}

// ************************************************************************* //

template<class CloudType>
void Foam::finiteRateNewtonLinODE<CloudType>::calculate
(
    const scalar t0, 
    const scalar dt,
    const scalar ha, 
    const scalar pc, 
    const scalar pT,
    const scalar pz,
    scalarField& Y
)
{
    if (pz < zLower_ || pz > zUpper_)
    {
//        cout << "return mf " << pz << "\n";
        return;
    }

//    cout << "calculate mf " << pz << "\n";
    label nspec = Y.size();
    
//    scalar spec0[nspec];
    scalar  *pSpec=&Y[0];
    
//    forAll(Y,np)
//    {
////        spec0[np] = Y[np];
////    cout << "scalar " << np << "\n";
////    cout << "Y= " << Y[np] << "\n";
////    cout << "spec= " << spec0[np] << "\n";
//    }
    
//    scalar Temp = pT;
    scalar Temp = pT;
//    cout << "Temp   " << Temp    << "\n";
    
    scalar pascalToCGS = 10.;
    scalar P_dynes = pc * pascalToCGS;
//    cout << "P_dynes" << P_dynes << "\n";
    
//    scalar delt = dt;
    scalar deltat = dt;
//    cout << "deltat " << deltat  << "\n";

    newtonlinodesolve_
    (
//        spec0,
        pSpec,
        &nspec,
        &Temp,
        &P_dynes,
        &deltat
    );

//    forAll(Y,np)
//    {
//        Y[np] = spec0[np];
//    cout << "scalar updated " << np << "\n";
//    cout << "spec= " << spec0[np] << "\n";
//    cout << "Y= " << Y[np] << "\n";
//    } 

}


// ************************************************************************* //

