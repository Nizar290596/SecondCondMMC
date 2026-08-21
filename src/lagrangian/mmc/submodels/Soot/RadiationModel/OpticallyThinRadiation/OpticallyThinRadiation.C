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

#include "OpticallyThinRadiation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::OpticallyThinRadiation<CloudType>::OpticallyThinRadiation
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudRadiationModel<CloudType>(dict,owner,typeName),
    TRef_(this->coeffDict().lookupOrDefault("TRef", 300.0))
{
    Info << "Environment Temperature for SP-Radiation Model = " << TRef_ << "K" << endl;
}


template<class CloudType>
Foam::OpticallyThinRadiation<CloudType>::OpticallyThinRadiation
(
    const OpticallyThinRadiation<CloudType>& cm
)
:
    CloudRadiationModel<CloudType>(cm),
    TRef_(this->coeffDict().lookupOrDefault("TRef", 300.0))
{
    Info << "Environment Temperature for SP-Radiation Model = " << TRef_ << "K" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::OpticallyThinRadiation<CloudType>::~OpticallyThinRadiation()
{}

// * * * * * * * * * * * * * * * * member function  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::OpticallyThinRadiation<CloudType>::Qrad
(
    const scalarField& Y,
    const scalar fVsoot,
    const scalar p,
    const scalar T
)
{
/*
    Model for emitting gas and soot species. No absorption.
    
    Gas species parameters set according to the baseline optically thin radiation model recommended
    by the TNF workshop http://www.sandia.gov/TNF/radiation.html
    
    Soot parameters set according to the baseline optically thin radiation model recommended by the 
    International Sooting Flames Workshop
    http://www.adelaide.edu.au/cet/isfworkshop/docs/isf_baseline_radiation_model.pdf
*/


    // Physical parameters
    const scalar environmentTemperature = TRef_;//300; // Kelvin
    
    const scalar Tenv4 = pow(environmentTemperature,4.);

    const scalar sigmaSB = 5.67e-8; // W/m^2K^4
    
    const scalar atm = 101.3e3; // Pa
    
    const scalar rho = this->owner().composition().particleMixture(Y).rho(p,T);


    // Find partial pressures in atmospheres; equals (mole fraction) * (pressure) / (standard atmospheric pressure)
    scalarField Pi(this->owner().composition().X(Y) * p / atm);
    
    // Find which radiative species are present in the mixture and assign their partial pressures
    scalar PH2O = 0;
    scalar PCO2 = 0;
    scalar PCO = 0;
    scalar PCH4 = 0;
    if(this->owner().composition().slgThermo().carrier().contains("H2O"))
        PH2O = Pi[this->owner().composition().slgThermo().carrier().species()["H2O"]];

    if(this->owner().composition().slgThermo().carrier().contains("CO2"))
        PCO2 = Pi[this->owner().composition().slgThermo().carrier().species()["CO2"]];

    if(this->owner().composition().slgThermo().carrier().contains("CO"))
        PCO = Pi[this->owner().composition().slgThermo().carrier().species()["CO"]];

    if(this->owner().composition().slgThermo().carrier().contains("CH4"))
        PCH4 = Pi[this->owner().composition().slgThermo().carrier().species()["CH4"]];

    // Absorption coefficients
    scalar uT = 1000./(T+SMALL);

    scalar K_H2O = -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));

    scalar K_CO2 =  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));

    scalar K_CO;
    if( T < 750. )  
        K_CO = 4.7869+T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)));
    else
        K_CO = 10.09+T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));

    // 4. Methane [1/m/bar]
    scalar K_CH4 = 6.6334 +T*(- 0.0035686+T*(1.6682e-08+T*(2.5611e-10-2.6558e-14*T)));

    scalar asTot = K_H2O*PH2O + K_CO2*PCO2 + K_CO*PCO + K_CH4*PCH4;
    
    // Soot
    scalar asSoot = 3.334e-4*fVsoot*T / (4.*sigmaSB);


    // Source term
    return  ( -4.*sigmaSB * (asTot * (pow(T,4.) - Tenv4) + asSoot * pow(T,4.)) ) / rho;     // Source term [W/m3]

}

// ************************************************************************* //

