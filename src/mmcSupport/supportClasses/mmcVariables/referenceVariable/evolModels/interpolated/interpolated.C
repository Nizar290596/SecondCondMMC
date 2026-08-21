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

#include "interpolated.H"

namespace Foam
{
    defineTypeNameAndDebug(interpolated, 0);
    addToRunTimeSelectionTable(evolModel,interpolated,word);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::interpolated::compute
(
    const particle& p,
    scalar dt, 
    const scalar& x
) const
{

    autoPtr<interpolation<scalar>> XiInterp
    (
        interpolation<scalar>::New
        (
            word
            (
                this->rVarProperties().subDict(refName_).lookup
                (
                    "interpolationScheme"
                )
            ),
            this->rVarProperties().db().lookupObject<volScalarField>(refName_)
        )
    );

    //return XiInterp().interpolate(p.position(),p.currentTetIndices());
    return XiInterp().interpolate(p.coordinates(),p.currentTetIndices());
}


inline Foam::scalar Foam::interpolated::initialise
(
    label patch,
    label patchFace,
    particle& p
) const
{
    return this->rVarProperties().db().lookupObject<volScalarField>
                        (refName_).boundaryField()[patch][patchFace];
}


inline Foam::scalar Foam::interpolated::initialise
(
    label celli,
    particle& p
) const
{
    return this->rVarProperties().db().lookupObject<volScalarField>
                        (refName_)[celli];
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolated::interpolated
(
    const word& refName,
    const fvMesh& mesh
)
:
    evolModel(mesh),

    refName_(refName)

//    XiRField_(mesh.objectRegistry::lookupObject<volScalarField>(refName))

//    XiInterp_
//    (
//        interpolation<scalar>::New
//        (
//            word
//            (
//                this->rVarProperties().subDict(refName).lookup
//                (
//                    "interpolationScheme"
//                )
//            ),
//            mesh.objectRegistry::lookupObject<volScalarField>(refName)
//        )
//    )
{}

Foam::interpolated::interpolated
(
    const interpolated& inpt
)
:
    evolModel(inpt),

    refName_(inpt.refName_)

//    XiInterp_
//    (
//        interpolation<scalar>::New
//        (
//            word
//            (
//                this->rVarProperties().subDict(refName_).lookup
//                (
//                    "interpolationScheme"
//                )
//            ),
//            inpt.XiInterp_().psi()
//        )
//    )
{}

// ************************************************************************* //
