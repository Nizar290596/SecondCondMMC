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

Class
    Foam::NCIrecycledTimeVaryingMappedFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.  Supplied data should be specified in
    constant/boundaryData/\<patchname\> where:
    - points : pointField with locations
    - ddd    : supplied values at time ddd
    The default mode of operation (mapMethod planarInterpolation) is
    to project the points onto a plane (constructed from the first threee
    points) and construct a 2D triangulation and finds for the face centres
    the triangle it is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and
    triangulation and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

Usage
    \table
        Property     | Description             | Required    | Default value
        setAverage   | flag to activate setting of average value | yes |
        perturb      | perturb points for regular geometries | no | 1e-5
        fieldTableName | alternative field name to sample | no| this field name
        mapMethod    | type of mapping | no | planarInterpolation
        offset   | for applying offset to mapped values  | no | constant 0.0
    \endtable

    \verbatim
    <patchName>
    {
        type            NCIrecycledTimeVaryingMappedFixedValue;
        setAverage      false;
        //perturb       0.0;
        //fieldTableName samples;
        //offset    constant 0.2;
        //mapMethod       nearest;
        setScale        true;
        ScaleFactorRMS  (1 1 4);
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    NCIrecycledTimeVaryingMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef NCIrecycledTimeVaryingMappedFixedValueFvPatchField_H
#define NCIrecycledTimeVaryingMappedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class NCIrecycledTimeVaryingMappedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NCIrecycledTimeVaryingMappedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories (folders)
        instantList sampleTimes_;

        //- Current starting index in sampleTimes (index of folders)
        label startSampleTime_;

        //- Interpolated values from startSampleTime 
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes (index of folders)
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Number of cycles done over boundaryData
        label nRecycles_;

        //- First instant stored in boundaryData
        scalar t0_;

        //- Time stored in boundaryData folder
        scalar DeltaT_;

        //- if true scale the rms of fields
        Switch setScale_;

        //- scaling factors
         const Type ScaleFactorRMS_;

        //- scaling factors
        const Type ScaleFactorMEAN_;

        //- mean field 
        Field<Type> meanSampledValues_;

        //start time value
        scalar subStartTimeValue_;

        //end time value
        scalar subEndTimeValue_;

        //Combined BC data file
        IOdictionary valsFile_;


public:

    //- Runtime type information
    TypeName("NCIrecycledTimeVaryingMappedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        //  onto a new patch
        NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        (
            const NCIrecycledTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        (
            const NCIrecycledTimeVaryingMappedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new NCIrecycledTimeVaryingMappedFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        NCIrecycledTimeVaryingMappedFixedValueFvPatchField
        (
            const NCIrecycledTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new NCIrecycledTimeVaryingMappedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return startSampledValues
            const Field<Type> startSampledValues()
            {
                 return startSampledValues_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Utility functions

            //- Find boundary data inbetween current time and interpolate
            void checkTable();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NCIrecycledTimeVaryingMappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
