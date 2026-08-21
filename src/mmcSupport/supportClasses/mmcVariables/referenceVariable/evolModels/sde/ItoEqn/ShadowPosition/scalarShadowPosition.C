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

#include "scalarShadowPosition.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Member Data * * * * * * * * * * * * * * //

#include "scalarShadowPositionModels.C"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ShadowPosition<Foam::scalar>::ShadowPosition
(
    const dictionary& entry,
    const fvMesh& mesh
)
:
    scalarItoEqn(),

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

    componentName_(entry.lookup("component")),

    component_
    (
        (componentName_ == "x") ? vector::X :
        (componentName_ == "y") ? vector::Y : vector::Z
    ),

    modelType_
    (
        word
        (
            IOdictionary
            (
                IOobject
                (
                    Foam::turbulenceModel::propertiesName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).lookup("simulationType")
        )
    ),

    intpSchemes_(entry.subDict("interpolationSchemes")),

//    UInterp_
//    (
//        interpolation<scalar>::New
//        (
//            intpSchemes_,
//            (mesh.objectRegistry::
//                lookupObject<volVectorField>("U")).component(component_)
//        )
//    ),

    UInterp_
    (
        interpolation<vector>::New
        (
            intpSchemes_,
            mesh.objectRegistry::lookupObject<volVectorField>("U")
        )
    ),

    DtInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            mesh.objectRegistry::lookupObject<volScalarField>("Dt")
        )
    ),

    DInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            mesh.objectRegistry::lookupObject<volScalarField>("D")
        )
    ),

    DeltaEInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            mesh.objectRegistry::lookupObject<volScalarField>("DeltaE")
        )
    ),
    
    vbInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            mesh.objectRegistry::lookupObject<volScalarField>("vb")
        )
    )
{
    constructTimeScaleModelTable();
}

Foam::ShadowPosition<Foam::scalar>::ShadowPosition
(
    const ShadowPosition& sp
)
:
    scalarItoEqn(),

    a_(sp.a_),

    b_(sp.b_),
    
    variableLambda_(sp.variableLambda_),

    componentName_(sp.componentName_),

    component_(sp.component_),

    modelType_(sp.modelType_),

    intpSchemes_(sp.intpSchemes_),

    UInterp_
    (
        interpolation<vector>::New
        (
            intpSchemes_,
            sp.UInterp_().psi()
        )
    ),

    DtInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            sp.DtInterp_().psi()
        )
    ),

    DInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            sp.DInterp_().psi()
        )
    ),

    DeltaEInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            sp.DeltaEInterp_().psi()
        )
    ),
    
    vbInterp_
    (
        interpolation<scalar>::New
        (
            intpSchemes_,
            sp.vbInterp_().psi()
        )         
    )
{}
// * * * * * * * * * * * * * *  Friend Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::LESomegaSP
(
    const vector& position,
    //const vector& coordinates,
    label cell,
    label face,
    const scalarShadowPosition& sp
)
{
    scalar Dt_p = sp.DtInterp_().interpolate(position,cell,face);
    //scalar Dt_p = sp.DtInterp_().interpolate(coordinates,cell,face);

    scalar D_p = sp.DInterp_().interpolate(position,cell,face);
    //scalar D_p = sp.DInterp_().interpolate(coordinates,cell,face);

    scalar DEff_p = Dt_p*0.4 + D_p*0.7;

    //- Filter size at particle location
    scalar DeltaE_p = sp.DeltaEInterp_().interpolate(position,cell,face);
    //scalar DeltaE_p = sp.DeltaEInterp_().interpolate(coordinates,cell,face);

    return (sp.a_* DEff_p / (sqr(DeltaE_p) + VSMALL));
}

Foam::scalar Foam::RASomegaSP
(
    const vector& position,
    //const vector& coordinates,
    label cell,
    label face,
    const scalarShadowPosition& sp
)
{
//    //- Effective diffusivity at particle location
//    scalar epsilon_p = epsilonInterp_().interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    //- Filter size at particle location
//    scalar k_p = kInterp_().interpolate
//    (
//        p.position(),p.cell(),p.face()
//    );

//    return ((a_/CT) * epsilon_p / k_p);
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::ShadowPosition<Foam::scalar>::constructTimeScaleModelTable()
{
    if (timeScaleModelTable_.size())
        return;

    timeScaleModelTable_.insert("LES",Foam::LESomegaSP);
    timeScaleModelTable_.insert("RAS",Foam::RASomegaSP);
}


Foam::scalar Foam::ShadowPosition<Foam::scalar>::omegaSPos
(
    const vector& position,
    //const vector& coordinates,
    label cell,
    label face
) const
{
    //typename HashTable<ShadowPosition<scalar>::timeScaleModelPtr>::iterator
    //    cstrIter = timeScaleModelTable_.find(modelType_);
    auto cstrIter = timeScaleModelTable_.find(modelType_);

   return cstrIter()(position,cell,face,*this);
   // return cstrIter()(coordinates,cell,face,*this);
}


void Foam::ShadowPosition<Foam::scalar>::computeDriftCoeff
(
    const particle& p,
    const scalar& x
) const
{
    //- Filtered velocity at particle location
    //scalar U_p = UInterp_().interpolate
    //                (p.position(),p.currentTetIndices())[component_];
    scalar U_p = UInterp_().interpolate
                    (p.coordinates(),p.currentTetIndices())[component_];

    //- relaxation frequency
    scalar omega = omegaSPos(p.position(),p.cell(),p.face());

    //scalar omega = omegaSPos(p.coordinates(),p.cell(),p.face());
    //- check if same from (Deff and Delta) or (epsilon/k)
//    scalar omegaSPos;
//    if (modelType_ == "LES")
//        omegaSPos = a_* DEff_p / sqr(DeltaE_p);
//    else if (modelType_ == "RANS")
//        omegaSPos = (a_/CT) * epsilon_p / k_p;

    scalar deltaT = UInterp_().psi().db().time().deltaTValue();

    scalar relaxT = 0;
    //- Relaxation of shadow towards particle location
    if (omega < 1e30)
    {
        scalar mixExtent = 1 - exp(-omega * deltaT);
        relaxT = (p.position()[component_]-x)*mixExtent;
	//	relaxT = (p.coordinates()[component_]-x)*mixExtent;
//        Info << "omega SP < 1e30 and relaxT = " << relaxT << endl;
        relaxT = relaxT/deltaT; 
    }

    //- Drift term of shadow Position
    this->A_ = U_p + relaxT;
}


void Foam::ShadowPosition<Foam::scalar>::computeDiffusionCoeff
(
    const particle& p
) const
{
    //- Effective diffusivity at particle location
//    scalar Dt_p = DtInterp_().interpolate
//    (
//        p.position(),p.currentTetIndices()
//    );

    //scalar D_p = DInterp_().interpolate
    //(
    //    p.position(),p.currentTetIndices()
    //);

    scalar D_p = DInterp_().interpolate
    (
        p.coordinates(),p.currentTetIndices()
    );

    scalar DEff_p = D_p ;//+ Dt_p;
    
    //scalar vb_p = vbInterp_().interpolate
    //(
    //    p.position(),p.currentTetIndices()
    //);

    scalar vb_p = vbInterp_().interpolate
    (
        p.coordinates(),p.currentTetIndices()
    );


    if (variableLambda_)
    {
	this->D_ = vb_p*sqrt(2.0*DEff_p);
	//Info <<"VariableLambda"<<endl;
    }
    else 
    {
    	this->D_ = b_*sqrt(2.0*DEff_p);
    }
}


Foam::scalar Foam::ShadowPosition<Foam::scalar>::t0
(
    const particle& p,
    scalar dt, 
    const scalar& x
) const
{
//    this->computeDiffusionCoeff(p);
//    scalar dw = rndGen_.Normal(0,1)*sqrt(dt);

    return p.position()[component_]; // + D_*dw
    //return p.coordinates()[component_];
}
// ************************************************************************* //
