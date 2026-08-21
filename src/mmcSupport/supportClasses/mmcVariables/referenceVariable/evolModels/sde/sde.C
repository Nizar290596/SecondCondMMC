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

#include "sde.H"

namespace Foam
{
    defineTypeNameAndDebug(sde, 0);
    addToRunTimeSelectionTable(evolModel,sde,word);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::sde::compute
(
//    TrackData& td,
    const particle& p,
    scalar dt, 
    const scalar& x
) const
{
    SDEModel_().computeDriftCoeff(p,x);

    SDEModel_().computeDiffusionCoeff(p);

    return SDEModel_().integrate(p,dt,x);
}


inline Foam::scalar Foam::sde::initialise
(
    label patch,
    label patchFace,
    particle& p
) const
{
    scalar dt = 1e-6;
    scalar x  = 0;
    return SDEModel_().t0(p,dt,x);
}


inline Foam::scalar Foam::sde::initialise
(
    label celli,
    particle& p
) const
{
    scalar dt = 1e-6;
    scalar x  = 0;
    return SDEModel_().t0(p,dt,x);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sde::sde
(
    const word& refName,
    const fvMesh& mesh
)
:
    evolModel(mesh),

    SDEModel_
    (
        ItoEqn<scalar>::New
        (
            this->rVarProperties().subDict(refName),
            mesh
        )
    )
{}


Foam::sde::sde
(
    const sde& sd
)
:
    evolModel(sd),

    SDEModel_(sd.SDEModel_().clone())
{}

// ************************************************************************* //
