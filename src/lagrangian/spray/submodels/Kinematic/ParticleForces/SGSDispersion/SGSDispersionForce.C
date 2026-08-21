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

\*---------------------------------------------------------------------------*/

#include "SGSDispersionForce.H"
#include "demandDrivenData.H"
#include "turbulenceModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SGSDispersionForce<CloudType>::CdRe(const scalar Re) const
{
    if (Re > 1000.0)
    {
        return 0.424*Re;
    }
    else
    {
        return 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));
    }
}

template<class CloudType>
Foam::vector Foam::SGSDispersionForce<CloudType>::dW(const scalar dt) const
{   
    Random& rnd = this->owner().rndGen();

    vector eta (rnd.GaussNormal<vector>());
    
    return eta*sqrt(dt);
}

template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::SGSDispersionForce<CloudType>::kModel() const
{
    const objectRegistry& obr = this->owner().mesh();
    const word turbName =
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            this->owner().U().group()
        );

    if (obr.foundObject<turbulenceModel>(turbName))
    {
        const turbulenceModel& model =
            obr.lookupObject<turbulenceModel>(turbName);
        return model.k();        
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SGSDispersionForce<CloudType>::SGSDispersionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    rndGen_(owner.rndGen()),
    C0_(readScalar(this->coeffs().lookup("C0"))),
    alpha_(this->coeffs().template lookupOrDefault<scalar>("alpha", 0.8)),
    tauTName_
    (
        this->coeffs().template lookupOrDefault<word>
        (
           "timeScale",
           "original"
        )
    ),
    Delta_ 
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>("delta")
    ),
    kPtr_(nullptr),
    ownK_(false)
{
    if (tauTName_ == "original" || tauTName_ == "velBased" || tauTName_ == "kSGSBased")
    {
	Info<< "    Choosing LES-SGS Disperson Model "
        << "with " << tauTName_ << " time scale model."
        << endl;    	
    }
    else
    {
        FatalErrorInFunction
        	<< " Invalid SGS Dispersion Time Scale Model : "
        	<< tauTName_ << nl 
            << "    Turblent Time scale must be either: " 
            << nl << "    'original', 'velBased' or 'kSGSBased'"
            << nl << exit(FatalError);
    }
}


template<class CloudType>
Foam::SGSDispersionForce<CloudType>::SGSDispersionForce
(
    const SGSDispersionForce& sgsdf
)
:
    ParticleForce<CloudType>(sgsdf),
    rndGen_(sgsdf.rndGen_),
    C0_(sgsdf.C0_),
    alpha_(sgsdf.alpha_),
    tauTName_(sgsdf.tauTName_),
    Delta_(sgsdf.Delta_),   
    kPtr_(nullptr),
    ownK_(false)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SGSDispersionForce<CloudType>::~SGSDispersionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SGSDispersionForce<CloudType>::cacheFields(const bool store)
{

    if (store)
    {
        tmp<volScalarField> tk = kModel();

        if (tk.isTmp())
        {
            kPtr_ = tk.ptr();
            ownK_ = true;
        }
        else
        {
            kPtr_ = &tk();
            ownK_ = false;
        }
    }
    else
    {

        if (ownK_ && kPtr_)
        {
            deleteDemandDrivenData(kPtr_);
            ownK_ = false;
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::SGSDispersionForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const label celli = p.cell();

    // Unresolved kinetic energy field [m2/s2]
    const volScalarField& k = *kPtr_;

    // Unresolved kinetic energy at particle location cell [m2/s2]
    const scalar kSGSc = k[celli];

    // Particle response time scale
    const scalar tauPInv = 0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()));
    const scalar tauP = 1.0/tauPInv;

    // LES Filter Width
    const scalar Delta = Delta_[celli];

    // Turbulent time scale
    // scalar tauT = 0;
    vector Source(Zero);   
    if (tauTName_ == "original")  // Original Bini and Jones (2007)
    {
/*
        scalar n = 2*alpha() - 1;
        if (n < 0)
        {
        	FatalErrorInFunction
                << "SGS dispersion model time scale exponent is -ve!" << nl
                << "Check if alpha is < 0.5." << nl
                << abort(FatalError);
        }
        // tauT = tauP * pow((tauP*sqrt(kSGSc)/Delta), n);
        Source = mass * sqrt(C0()*kSGSc) * dW(dt) * sqrt(pow(Delta,n)) / (sqrt(tauP*pow(tauP*sqrt(kSGSc),n))*dt + ROOTVSMALL);     	
*/      
	scalar tauT = pow(tauP,2*alpha())/(pow(Delta/sqrt(kSGSc),2*alpha()-1));   
        scalar b = sqrt(C0()*kSGSc/tauT);
        Source = mass/dt * b * dW(dt);
    }
    else if (tauTName_ == "velBased") // Based on particle velocity
	{
    	// tauT = Delta / (mag(p.U()) + ROOTVSMALL);
        Source = mass * sqrt(C0()*kSGSc) * dW(dt) * sqrt(mag(p.U())) / (sqrt(Delta)*dt + ROOTVSMALL);

    }
    else if (tauTName_ == "kSGSBased") // Based on unresolved K
	{
    	// tauT =  Delta / (sqrt(kSGSc) + ROOTVSMALL);
        Source = mass * sqrt(C0()*kSGSc) * dW(dt) * sqrt(sqrt(kSGSc)) / (sqrt(Delta)*dt + ROOTVSMALL);    	
    }
  	
    // Explicit Non-Coupled Force
    // value.Su() = mass*sqrt(C0()*kSGSc/(tauT + VSMALL))*dW(dt)/dt;
    value.Su() = Source;

    return value;
}

// ************************************************************************* //
