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

Description

    Inlet with synthetic turbulent fluctuations (version 2).

    Uses the syntheticTurbulence class to generate random fluctuations. 
    The algorithm is contained in Foam::syntheticTurbulence.
    
    Integral length scale is specified via Foam::turbulenceLengthScale class
    and its derived.

    Example of use (file 0/U)
    
    \verbatim
    
    inlet
    {
        type    syntheticTurbulenceInlet;
    
        // parameters for synthetic turbulence generation
        // (see also Foam::syntheticTurbulence class)

        kName            k;         // TKE
        omegaName        omega;     // specific dissipaton rate of TKE
        UrefName         Uref;      // reference (average) velocity

        distType         uniform;   // uniform (default)|logarithmic

        betaStar         0.09;      // epsilon = betaStar*k*omega

        A                1.4526;

        ft               1.0;    // time scale factor (1 if not specified)
        fampl            1.0;    // amplitude factor (1 if not specified)
        fk               1.0;    // min. wavenumber factor (1 if not specified)
        fSpect           1.0;    // spectrum multiplier (1 if not specified)
        nModes           100;    // number of frequency modes to generate
        minK             0.0;    // min K (turb. kin. energy) to switch on S.T.
        kwMin            0.0;    // min wavenumber
        kwMax            1.0;    // max wavenumber
                                 // if not specified the max cell size is 
                                 // used for calculation

        dtRatio          1;      // ratio between computational dt and
                                 // timestep used for temporal filter 
                                 // (1 if not specified)

        //specification of integral length scale (see related class)
        L
        {
            type        kwLengthScale;
            kName       k;
            omegaName   omega;
            fL          1;
        }

        // wall blending function parameters
        bScaling        1;      // slope factor of blending function
        minfbl          0.1;    // limiting value of blending function
        deltaIn         0.126;  //boundary layer thickness at inlet

        // placeholder
        value           uniform ( 0 0 0 );

    }
    
    
    
    \endverbatim

SourceFiles
    syntheticTurbulenceInletFvPatchField.C


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

#ifndef syntheticTurbulenceInletFvPatchField_H
#define syntheticTurbulenceInletFvPatchField_H

//v2
#include "syntheticTurbulence.H"

#include "Pstream.H"
#include "fvCFD.H"
#include "Random.H"
#include "fixedValueFvPatchField.H"
#include "volMesh.H"
#include "volFields.H"
#include "Switch.H"
//#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class syntheticTurbulenceInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class syntheticTurbulenceInletFvPatchField
:
    public fixedValueFvPatchField<vector>
{
    // Private data

        syntheticTurbulence* turbulencePtr_;

        // fields
        scalarField y_;
        scalarField fbl_;
        
        // parameters for wall blending function
        scalar deltaIn_;
        scalar bScaling_;
        scalar minfbl_;
        
    // private functions

public:

    //- Runtime type information
    TypeName("syntheticTurbulenceInletUMeanControl");

    virtual ~syntheticTurbulenceInletFvPatchField();

    // Constructors

        //- Construct from patch and internal field
        syntheticTurbulenceInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        syntheticTurbulenceInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //    syntheticTurbulenceInletFvPatchField onto a new patch
        syntheticTurbulenceInletFvPatchField
        (
            const syntheticTurbulenceInletFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        syntheticTurbulenceInletFvPatchField
        (
            const syntheticTurbulenceInletFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new syntheticTurbulenceInletFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        syntheticTurbulenceInletFvPatchField
        (
            const syntheticTurbulenceInletFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new syntheticTurbulenceInletFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
