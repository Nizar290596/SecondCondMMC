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
    Foam::synteticTurbulence


Authors
    Andrea Montorfano, Federico Piscaglia
    Dipartimento di Energia, Politecnico di Milano
    via Lambruschini 4
    I-20156 Milano (MI)
    ITALY

Contact
    andrea.montorfano@polimi.it , ph. +39 02 2399 8639
    federico.piscaglia@polimi.it, ph. +39 02 2399 8620    


        
\*---------------------------------------------------------------------------*/

// Modified von Karman energy spectrum. Cfr. eq. (4) on:
// L. Davidson, "using isotropic synthetic fluctuations as inlet 
// boundary conditions for unsteady simulations", 
// Advances and Applications in Fluid Mechanics, 2007
// Note: 3*uRms^2 = 2*k_

inline scalarField syntheticTurbulence::MVKS_
(
    const label i,
    const scalar ke,
    const scalarField kw,
    const scalar kkol
)
{
    
    scalar c1 = A_/ke*2.0*k_[i]/3;

    scalarField den = pow((1.0+sqr(kw/ke)),(17.0/6.0));

    scalarField c2 =  pow((kw/ke),4.0)/den;

    scalarField c3 =  exp(-2.0*sqr(kw/kkol));

    scalarField Ek = c1*c2*c3;
    
    return  Ek;
}


inline void syntheticTurbulence::logDist_
(
    scalarField& kw
)
{
    if (kw[0] > SMALL && kwMax_ > SMALL)
    { 
        scalar deltakl = (log(kwMax_) - log(kw[0])) / (nModes_-1);
        forAll(kw,n)
        {
            kw[n] = exp(log(kw[0])+ n*deltakl);
        }
    }    
    else
    {
        forAll(kw,n)
        {
            kw[n] = 0.0;
        }
    }
    return;
}
    

inline void syntheticTurbulence::uniDist_
(
    scalarField& kw
)
{
    if (kw[0] > SMALL && kwMax_ > SMALL)
    { 
        scalar deltakl = (kwMax_ - kw[0]) / (nModes_-1);
        forAll(kw,n)
        {
            kw[n] =  kw[0]+n*deltakl;
        }
    }    
    else
    {
        forAll(kw,n)
        {
            kw[n] = 0.0;
        }
    }
}


inline Foam::vectorField Foam::syntheticTurbulence::blend
(
    vectorField& ut
)
{
    return ut*fbl_;
}


template<class Type>
Foam::Field<Type> Foam::syntheticTurbulence::readPatchField
(
    const word& fieldName, 
    const fvPatchVectorField& pvf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    Field<Type> pf(coords_.size(),pTraits<Type>::zero);
    
    if(debug)
    {
        Info << "retrieving field: " << fieldName << endl;
    }
    
    if 
    (
        ! runTime_.foundObject<vfType>(fieldName)
    )
    {

        if(debug)
        {
            Info << fieldName 
                << " is not in the database. Read from file" << endl;
        }

         regIOobject::store
        (
            new vfType
            (
                IOobject
                (
                    fieldName,
                    runTime_.timeName(),
                    runTime_,    
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        
    }
    else
    {
        if (debug)
        {
            Info << fieldName
                << " is in the database. Load from memory" << endl;
        }
    }
    
    const vfType& fieldRef
    (
        runTime_.lookupObject<vfType>(fieldName)
    );

    pf = fieldRef.boundaryField()[pvf.patch().index()];

    return pf;

}
