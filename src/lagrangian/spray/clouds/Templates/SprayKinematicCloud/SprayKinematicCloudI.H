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

#include "fvmSup.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::SprayKinematicCloud<CloudType>&
Foam::SprayKinematicCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::SprayKinematicCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::SprayKinematicCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::SprayKinematicCloud<CloudType>::fuelCloud1Properties() const
{
    return fuelCloud1Properties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::SprayKinematicCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::SprayKinematicCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::SprayKinematicCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::SprayKinematicCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::SprayKinematicCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::SprayKinematicCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::SprayKinematicCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::SprayKinematicCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::SprayKinematicCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::SprayKinematicCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::SprayKinematicCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::SprayKinematicCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::SprayKinematicCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
inline const typename Foam::SprayKinematicCloud<CloudType>::forceType&
Foam::SprayKinematicCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::SprayKinematicCloud<CloudType>::forceType&
Foam::SprayKinematicCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::SprayKinematicCloud<CloudType>::functionType&
Foam::SprayKinematicCloud<CloudType>::functions()
{
    return functions_;
}


// Access to sub models 
template<class CloudType>
inline const Foam::InjectionModelList<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::injectors() const
{
    return injectors_;
}
 
 
template<class CloudType>
inline Foam::InjectionModelList<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::injectors()
{
    return injectors_;
}
 
 
template<class CloudType>
inline const Foam::DispersionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::dispersion() const
{
    return *dispersionModel_;
}
 
 
template<class CloudType>
inline Foam::DispersionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::dispersion()
{
    return *dispersionModel_;
}
 
 
template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::patchInteraction() const
{
    return *patchInteractionModel_;
}
 
 
template<class CloudType>
inline Foam::PatchInteractionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::patchInteraction()
{
    return *patchInteractionModel_;
}
 
 
template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::stochasticCollision() const
{
    return *stochasticCollisionModel_;
}
 
 
template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::stochasticCollision()
{
    return *stochasticCollisionModel_;
}
 
 
template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::surfaceFilm() const
{
    return *surfaceFilmModel_;
}
 
 
template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::surfaceFilm()
{
    return *surfaceFilmModel_;
}
 
 
template<class CloudType>
inline const Foam::PackingModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::packingModel() const
{
    return *packingModel_;
}
 
 
template<class CloudType>
inline Foam::PackingModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::packingModel()
{
    return *packingModel_;
}
 
 
template<class CloudType>
inline const Foam::DampingModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::dampingModel() const
{
    return *dampingModel_;
}
 
 
template<class CloudType>
inline Foam::DampingModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::dampingModel()
{
    return *dampingModel_;
}
 
 
template<class CloudType>
inline const Foam::IsotropyModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::isotropyModel() const
{
    return *isotropyModel_;
}
 
 
template<class CloudType>
inline Foam::IsotropyModel<Foam::SprayKinematicCloud<CloudType>>&
Foam::SprayKinematicCloud<CloudType>::isotropyModel()
{
    return *isotropyModel_;
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::SprayKinematicCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::label Foam::SprayKinematicCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::SprayKinematicCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIters(*this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::SprayKinematicCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::SprayKinematicCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::SprayKinematicCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::SprayKinematicCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


template<class CloudType>
inline Foam::Random& Foam::SprayKinematicCloud<CloudType>::rndGen() const
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*>>&
Foam::SprayKinematicCloud<CloudType>::cellOccupancy()
{
    if (!cellOccupancyPtr_.good())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::SprayKinematicCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::SprayKinematicCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::SprayKinematicCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::SprayKinematicCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::SprayKinematicCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::scalar Foam::SprayKinematicCloud<CloudType>::dMin() const
{
    return dMin_;
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::SprayKinematicCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const volScalarField::Internal
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::SprayKinematicCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "vDotSweep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep.ref();
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        vDotSweep[celli] += p.nParticle()*p.areaP()*mag(p.U() - U_[celli]);
    }

    vDotSweep.primitiveFieldRef() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::SprayKinematicCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta.ref();
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        theta[celli] += p.nParticle()*p.volume();
    }

    theta.primitiveFieldRef() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::SprayKinematicCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha.ref().primitiveFieldRef();
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        alpha[celli] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::SprayKinematicCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff.ref().primitiveFieldRef();
    forAllConstIters(*this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        rhoEff[celli] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


template<class CloudType>
inline Foam::eulerianStatistics& Foam::SprayKinematicCloud<CloudType>::eulerianStats()
{
    return eulerianStats_;
}


template<class CloudType>
inline const Foam::eulerianStatistics& Foam::SprayKinematicCloud<CloudType>::eulerianStats() const
{
    return eulerianStats_;
}


template<class CloudType>
inline const Foam::dictionary& Foam::SprayKinematicCloud<CloudType>::eulerianStatsDict() const
{
    return eulerianStatsDict_;
}


// ************************************************************************* //
