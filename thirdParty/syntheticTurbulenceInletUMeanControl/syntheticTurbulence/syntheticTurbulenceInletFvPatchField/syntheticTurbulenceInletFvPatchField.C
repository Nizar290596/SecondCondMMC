/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This file is not part of OpenFOAM but an original routine developed 
    in the OpenFOAM techonology.

    You can redistribute it and/or modify it under the terms of the GNU General
    Public License as published by the Free Software Foundation; either version
    3 of the License, or (at your option) any later version. See GNU General
    Public License at <http://www.gnu.org/licenses/gpl.html>

    OpenFOAM is a trademark of OpenCFD.

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


Class
    Foam::syntheticTurbulenceInletFvPatchField

Authors
    Andrea Montorfano, Federico Piscaglia
    Dipartimento di Energia, Politecnico di Milano
    via Lambruschini 4
    I-20156 Milano (MI)
    ITALY

Contact
    andrea.montorfano@polimi.it , ph. +39 02 2399 3909
    federico.piscaglia@polimi.it, ph. +39 02 2399 8620    

    
            
\*---------------------------------------------------------------------------*/

#include "syntheticTurbulenceInletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

syntheticTurbulenceInletFvPatchField::syntheticTurbulenceInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    turbulencePtr_(NULL),
    y_(p.size()),
    fbl_(p.size()),
    deltaIn_(0.0),
    bScaling_(0.0),
    minfbl_(0.0)
{}


syntheticTurbulenceInletFvPatchField::syntheticTurbulenceInletFvPatchField
(
    const syntheticTurbulenceInletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    turbulencePtr_(ptf.turbulencePtr_),
    y_(ptf.y_),
    fbl_(ptf.fbl_),
    deltaIn_(ptf.deltaIn_),
    bScaling_(ptf.bScaling_),
    minfbl_(ptf.minfbl_)
{
    turbulencePtr_->operator++();
}

// ----------------------------------------------------------------------
// explicit constructor

syntheticTurbulenceInletFvPatchField::syntheticTurbulenceInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    turbulencePtr_(NULL),
    y_(p.size()),
    fbl_(p.size()),
    deltaIn_(readScalar(dict.lookup("deltaIn"))),
    bScaling_(readScalar(dict.lookup("bScaling"))),
       minfbl_(readScalar(dict.lookup("minfbl")))
{

    fixedValueFvPatchField<vector>::operator==
    (
        Field<vector>("value", dict, p.size())
    );

// retrieve mesh reference and read utf (if present)
    const fvMesh& mesh = patch().boundaryMesh().mesh(); 

    // debug
    if (!db().foundObject<volVectorField>("utf"))
    {
         regIOobject::store
        (
            new volVectorField
            (
                IOobject
                (
                    "utf",
                    mesh.time().timeName(),
                    this->db(),    
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("utf",dimVelocity,pTraits<vector>::zero)
            )
        );
    }

    const volVectorField& 
        volUtf=db().lookupObject<volVectorField>("utf");

    // calculate scalarFields
    //wallDist y(mesh);
        //volScalarField y(wallDist::New(mesh).y());
        // volScalarField& y(wallDist::New(mesh));
           //y_ = y.boundaryField()[this->patch().index()];

    fbl_ = 1.0; // DAMPING TURNED OFF!!
    //fbl_ = max
    //(
    //    0.5*(1.0-Foam::tanh((y_ - deltaIn_)/bScaling_)),
    //    minfbl_
    //);        
    
    

    // create syntheticTurbulence
    turbulencePtr_ = new syntheticTurbulence 
    (
        dict,
        *this,
        volUtf.boundaryField()[this->patch().index()],
        fbl_
    );

}

// ----------------------------------------------------------------------

syntheticTurbulenceInletFvPatchField::syntheticTurbulenceInletFvPatchField
(
    const syntheticTurbulenceInletFvPatchField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    turbulencePtr_(ptf.turbulencePtr_),
    y_(ptf.y_),
    fbl_(ptf.fbl_),
    deltaIn_(ptf.deltaIn_),
    bScaling_(ptf.bScaling_),    
    minfbl_(ptf.minfbl_)
{
    turbulencePtr_->operator++();
}


syntheticTurbulenceInletFvPatchField::syntheticTurbulenceInletFvPatchField
(
    const syntheticTurbulenceInletFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    turbulencePtr_(ptf.turbulencePtr_),
    y_(ptf.y_),
    fbl_(ptf.fbl_),
    deltaIn_(ptf.deltaIn_),
    bScaling_(ptf.bScaling_),    
    minfbl_(ptf.minfbl_)
{
    turbulencePtr_->operator++();
}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void syntheticTurbulenceInletFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
}


void syntheticTurbulenceInletFvPatchField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);
}


void syntheticTurbulenceInletFvPatchField::write(Ostream& os) const
{

    fvPatchField<vector>::write(os);

    os.writeKeyword("minfbl")       << minfbl_ << token::END_STATEMENT << nl;
    os.writeKeyword("deltaIn")       << minfbl_ << token::END_STATEMENT << nl;
    os.writeKeyword("bScaling")       << minfbl_ << token::END_STATEMENT << nl;        

    // write syntheticTurbulence dict
    turbulencePtr_->write(os);    
    
    this->writeEntry("value", os);

}

syntheticTurbulenceInletFvPatchField::~syntheticTurbulenceInletFvPatchField()
{
    //if (turbulencePtr_->okToDelete())
        if (turbulencePtr_->unique())
    {
        delete turbulencePtr_;
    }
    else
    {
        turbulencePtr_->operator--();
    }
}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::syntheticTurbulenceInletFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    turbulencePtr_->update();
    
    operator==(turbulencePtr_->U());

    fixedValueFvPatchField<vector>::updateCoeffs();

}

 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    syntheticTurbulenceInletFvPatchField
);

} // End namespace Foam

// ************************************************************************* //
