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

#include "physicoChemicalConstants.H"

using namespace Foam; //::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::SprayThermoCloud<CloudType>&
Foam::SprayThermoCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::SprayThermoCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::SprayThermoCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::SLGThermo& Foam::SprayThermoCloud<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::SprayThermoCloud<CloudType>::T() const
{
    return T_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::SprayThermoCloud<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::EvaporationModel<Foam::SprayThermoCloud<CloudType>>&
Foam::SprayThermoCloud<CloudType>::evaporation() const
{
    return evaporationModel_;
}


template<class CloudType>
inline Foam::EnvelopeModel<Foam::SprayThermoCloud<CloudType>>&
Foam::SprayThermoCloud<CloudType>::envelope()
{
    return envelopeModel_();
}

template<class CloudType>
inline const Foam::CpModel<Foam::SprayThermoCloud<CloudType> >& 
Foam::SprayThermoCloud<CloudType>::cpModel() const
{
    return cpModel_;
}

template<class CloudType>
inline const Foam::FPRadiationModel<Foam::SprayThermoCloud<CloudType> >& 
Foam::SprayThermoCloud<CloudType>::fpRadiationModel() const
{
    return fpRadiationModel_;
}

template<class CloudType>
inline const Foam::integrationScheme&
Foam::SprayThermoCloud<CloudType>::TIntegrator() const
{
    return TIntegrator_;
}

template<class CloudType>
inline const Foam::integrationScheme&
Foam::SprayThermoCloud<CloudType>::MIntegrator() const
{
    return MIntegrator_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::SprayThermoCloud<CloudType>::f() const
{
    return f_;
}

 
template<class CloudType>
inline const Foam::volScalarField& Foam::SprayThermoCloud<CloudType>::magSqrGradf() const
{
    return magSqrGradf_;
}

template<class CloudType>
inline Foam::StochasticLib1& Foam::SprayThermoCloud<CloudType>::rand_f_Gen()
{
    return rand_f_Gen_;
}

//- Return const access to the popeCloud
template<class CloudType>
inline Foam::basicAerosolReactingPopeCloud& Foam::SprayThermoCloud<CloudType>::pope() const
{
    return pope_;
}


template<class CloudType>
inline Foam:: List<label> Foam::SprayThermoCloud<CloudType>::specieId() const
{
    return specieId_;
}


template<class CloudType>
inline Foam:: label Foam::SprayThermoCloud<CloudType>::nFuelSpecies() const
{
    return nFuelSpecies_;
}


template<class CloudType>
inline Foam:: scalarField Foam::SprayThermoCloud<CloudType>::linkFG() const
{
    return linkFG_;
}


template<class CloudType>
inline Foam:: label Foam::SprayThermoCloud<CloudType>::nProductSpecies() const
{
    return nProductSpecies_;
}


template<class CloudType>
inline Foam:: scalarField Foam::SprayThermoCloud<CloudType>::linkPG() const
{
    return linkPG_;
}

template<class CloudType>
inline Foam:: scalarField Foam::SprayThermoCloud<CloudType>::YProducts() const
{
    return YProducts_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::SprayThermoCloud<CloudType>::MTrans()
{
    return MTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::SprayThermoCloud<CloudType>::MTrans() const
{
    return MTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::SprayThermoCloud<CloudType>::MTransRate()
{
    return MTransRate_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::SprayThermoCloud<CloudType>::MTransRate() const
{
    return MTransRate_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::SprayThermoCloud<CloudType>::SuMass(volVectorField& U) const
{
    if (debug)
    {
        Info<< "MTrans min/max = " << min(MTrans()).value() << ", "
            << max(MTrans()).value() << endl;
            
    }

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("U"))
        {
            notImplemented
            (
                "SemiImlicit code is not implemented"                
            );
             
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -MTrans()/ (this->db().time().deltaT());           

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::SprayThermoCloud<CloudType>::SM(volScalarField& rho) const
{
    if (debug)
    {
        Info<< "MTransRate min/max = " << min(MTransRate()).value() << ", "
            << max(MTransRate()).value() << nl;
    }

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("rho"))
        {
            notImplemented
            (
                "SemiImlicit code is not implemented"                
            );
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -MTransRate()/(this->db().time().deltaTValue());

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(rho, dimMass/dimTime));
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::SprayThermoCloud<CloudType>::SM() const
{
   tmp<volScalarField::Internal> trhoTrans
   (
       new volScalarField::Internal
       (
           IOobject
           (
               this->name() + ":rhoTrans",
               this->db().time().timeName(),
               this->db(),
               IOobject::NO_READ,
               IOobject::NO_WRITE,
               false
           ),
           this->mesh(),
           dimensionedScalar
           (
               "zero",
               dimMass/dimTime/dimVolume,
               0.0
           )
       )
   );    
   
   if (this->solution().coupled())
   {
       scalarField& sourceField = trhoTrans.ref();
       // forAll(MTransRate(), i)
       // {
       //     sourceField += MTransRate();
       // } 
       sourceField = MTransRate()/(this->db().time().deltaTValue()*this->mesh().V());
   } 

   return trhoTrans;
}


template<class CloudType>
inline Foam::List<int> Foam::SprayThermoCloud<CloudType>::parcelId() const 
{
    return parcelId_;
}

template<class CloudType>
inline Foam::List<int>& Foam::SprayThermoCloud<CloudType>::parcelId() 
{
    return parcelId_;
}


template<class CloudType>
inline Foam::scalar Foam::SprayThermoCloud<CloudType>::Tmax() const
{
    scalar T = -GREAT;
    scalar n = 0;
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        T = max(T, p.T());
        n++;
    }

    reduce(T, maxOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
inline Foam::scalar Foam::SprayThermoCloud<CloudType>::Tmin() const
{
    scalar T = GREAT;
    scalar n = 0;
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        T = min(T, p.T());
        n++;
    }

    reduce(T, minOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
inline Foam::List<int> Foam::SprayThermoCloud<CloudType>::emptySuperCell() const 
{
    return emptySuperCell_;
}

template<class CloudType>
inline Foam::List<int>& Foam::SprayThermoCloud<CloudType>::emptySuperCell() 
{
   
    return emptySuperCell_;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::SprayThermoCloud<CloudType>::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );
 
 
    return tEp;
}
 
 
template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::SprayThermoCloud<CloudType>::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );
 
 
    return tap;
}
 
 
template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::SprayThermoCloud<CloudType>::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );
 
 
    return tsigmap;
}


// ************************************************************************* //
