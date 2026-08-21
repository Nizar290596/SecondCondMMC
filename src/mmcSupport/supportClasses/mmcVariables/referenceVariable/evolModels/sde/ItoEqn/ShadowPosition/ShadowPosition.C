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

#include "ShadowPosition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class T>
Foam::ShadowPosition<T>::ShadowPosition
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    ItoEqn<T>(),//entry,mesh)

    a_
    (
        readScalar(entry.lookup("a"))
    ),

    b_
    (
        readScalar(entry.lookup("b"))
    ),
    
    variableLambda_
    (
        readBool(entry.lookup("variableLambda"))
    ),

    UInterp_
    (
        interpolation<vector>::New
        (
            word(entry.lookup("interpolationScheme")),
            mesh.objectRegistry::
                lookupObject<GeometricField<T, fvPatchField, volMesh>>("U")
        )
    ),

    DEffInterp_
    (
        interpolation<scalar>::New
        (
            word(entry.lookup("interpolationScheme")),
            mesh.objectRegistry::lookupObject<volScalarField>("DEff")
        )
    ),

    DeltaEInterp_
    (
        interpolation<scalar>::New
        (
            word(entry.lookup("interpolationScheme")),
            mesh.objectRegistry::lookupObject<volScalarField>("DeltaE")
        )
    ),
    
    vbInterp_
    (
        interpolation<scalar>::New
        (
            word(entry.lookup("interpolationScheme")),
            mesh.objectRegistry::lookupObject<volScalarField>("vb")
        )
    )
    
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::ShadowPosition<T>::computeDriftCoeff
(
    const T& x,
    const particle& p
)
{

//    NOT IMPLEMENTED!!!!!!!!!!!!!!!!
//    //- Filtered velocity at particle location
//    vector U_p = UInterp_.interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    //- Effective diffusivity at particle location
//    scalar DEff_p = DEffInterp_.interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    //- Filter size at particle location
//    scalar DeltaE_p = DeltaEInterp_.interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    //- relaxation frequency
//    scalar wSPos = a_* DEff_p / sqr(DeltaE_p);

//    //- Relaxation of shadow towards particle location
//    vector relaxT = (x - p.position())*wSPos;

//    //- Drift term of shadow Position
//    this->A() = U_p + relaxT;
}

template<class T>
inline void Foam::ShadowPosition<T>::computeDiffusionCoeff(const particle& p)
{
//    NOT IMPLEMENTED!!!!!!!!!!!!!!!!
//    //- Effective diffusivity at particle location
//    scalar DEff_p = DEffInterp_.interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    this->D() = sqrt(2.0*DEff_p);
}
// ************************************************************************* //
