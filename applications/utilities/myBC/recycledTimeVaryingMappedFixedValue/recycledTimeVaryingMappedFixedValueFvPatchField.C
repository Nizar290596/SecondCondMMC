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

#include "recycledTimeVaryingMappedFixedValueFvPatchField.H"
#include "Time.H"
#include "AverageField.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::
recycledTimeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    setAverage_(false),
    perturb_(0),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(),
    nRecycles_(0),
    t0_(0),
    DeltaT_(0),
    setScale_(0),
    ScaleFactorRMS_(Zero),
    ScaleFactorMEAN_(1.0),
    meanSampledValues_(0)
{}


template<class Type>
Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::
recycledTimeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    fieldTableName_(iF.name()),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
//    offset_(Function1<Type>::New("offset", dict)),
    nRecycles_(0),
    t0_(0),
    DeltaT_(0),
    setScale_(dict.lookupOrDefault("setScale", false)),
    ScaleFactorRMS_(dict.lookupOrDefault<Type>("ScaleFactorRMS", Zero)),
    ScaleFactorMEAN_(dict.lookupOrDefault<scalar>("ScaleFactorMEAN", 1.0)),
    meanSampledValues_(0)
{
    Info << " scalingFactorRMS = " << ScaleFactorRMS_<<endl;
    Info << " scalingFactorMEAN = " << ScaleFactorMEAN_<<endl;

    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }


    dict.readIfPresent("fieldTable", fieldTableName_);

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
//        this->evaluate(Pstream::blocking);
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::
recycledTimeVaryingMappedFixedValueFvPatchField
(
    const recycledTimeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(ptf.offset_, false),
    nRecycles_(ptf.nRecycles_),
    t0_(ptf.t0_),
    DeltaT_(ptf.DeltaT_),
    setScale_(ptf.setScale_),
    ScaleFactorRMS_(ptf.ScaleFactorRMS_),
    ScaleFactorMEAN_(ptf.ScaleFactorMEAN_),
    meanSampledValues_(0)
{}


template<class Type>
Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::
recycledTimeVaryingMappedFixedValueFvPatchField
(
    const recycledTimeVaryingMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false),
    nRecycles_(ptf.nRecycles_),
    t0_(ptf.t0_),
    DeltaT_(ptf.DeltaT_),
    setScale_(ptf.setScale_),
    ScaleFactorRMS_(ptf.ScaleFactorRMS_),
    ScaleFactorMEAN_(ptf.ScaleFactorMEAN_),
    meanSampledValues_(ptf.meanSampledValues_)
{}


template<class Type>
Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::
recycledTimeVaryingMappedFixedValueFvPatchField
(
    const recycledTimeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false),
    nRecycles_(ptf.nRecycles_),
    t0_(ptf.t0_),
    DeltaT_(ptf.DeltaT_),
    setScale_(ptf.setScale_),
    ScaleFactorRMS_(ptf.ScaleFactorRMS_),
    ScaleFactorMEAN_(ptf.ScaleFactorMEAN_),
    meanSampledValues_(ptf.meanSampledValues_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    if (startSampledValues_.size())
    {
        startSampledValues_.autoMap(m);
        endSampledValues_.autoMap(m);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const recycledTimeVaryingMappedFixedValueFvPatchField<Type>& tiptf =
        refCast<const recycledTimeVaryingMappedFixedValueFvPatchField<Type>>(ptf);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::checkTable()
{
    // Initialise
    if (mapperPtr_.empty())
    {
//        pointIOField samplePoints
//        (
//            IOobject
//            (
//                "points",
//                this->db().time().constant(),
//                "boundaryData"/this->patch().name(),
//                this->db(),
//                IOobject::MUST_READ,
//                IOobject::AUTO_WRITE,
//                false
//            )
//        );

//        const fileName samplePointsFile = samplePoints.filePath();

        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().constant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

        if (debug)
        {
            Info<< "recycledTimeVaryingMappedFixedValueFvPatchField :"
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
        }

        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                 this->patch().patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );

        // Read the times for which data is available
        const fileName samplePointsDir = samplePointsFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);

        // Set first time stored in boundarayData & time stored
        t0_ = sampleTimes_.first().value();
        DeltaT_ = sampleTimes_.last().value() - t0_;

        if (debug)
        {
            Info<< "recycledTimeVaryingMappedFixedValueFvPatchField : In directory "
                << samplePointsDir << " found times "
                << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
                << endl;
        }

        if (setScale_)
        {
             fileName sampleMeanFile
            (
                this->db().time().constant()
               /"boundaryData"
               /this->patch().name()
               /fieldTableName_+"Mean"
            );

            //Field<Type> MeanField((IFstream(sampleMeanFile)()));

            Info<< "recycledTimeVaryingMappedFixedValueFvPatchField :"
                << " Read sample mean from "
                << sampleMeanFile << endl;

            Field<Type> MeanField;
            //Info << "mean field element 1 : " << MeanField[0] << endl;

            IFstream(sampleMeanFile)() >> MeanField;

            if (MeanField.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << MeanField.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << sampleMeanFile << exit(FatalError);
            }

            meanSampledValues_ = mapperPtr_().interpolate(MeanField);
            Info << "size of interpolated meanfield = " << meanSampledValues_.size()<< endl;

        }

        
    }


    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    scalar timeVal = this->db().time().value();
    label n = 0;

    if (DeltaT_ > 0)
    {
        n = label((timeVal-t0_)/DeltaT_);
        timeVal = timeVal - n*DeltaT_;
    }

    if (n > nRecycles_)
    {
        startSampleTime_ = -1;
        nRecycles_ = n;
    }

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        timeVal,
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << this->db().time().value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }

    Info << this->patch().name() << " lowTime = " << sampleTimes_[lo].name()
            << " highTime = " << sampleTimes_[hi].name() << endl;

    scalar tstart = this->db().time().elapsedCpuTime();

    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[lo].name()
                    << endl;
            }


//            // Reread values and interpolate
//            AverageIOField<Type> vals
//            (
//                IOobject
//                (
//                    fieldTableName_,
//                    this->db().time().constant(),
//                    "boundaryData"
//                   /this->patch().name()
//                   /sampleTimes_[startSampleTime_].name(),
//                    this->db(),
//                    IOobject::MUST_READ,
//                    IOobject::AUTO_WRITE,
//                    false
//                )
//            );

//            if (vals.size() != mapperPtr_().sourceSize())
//            {
//                FatalErrorInFunction
//                    << "Number of values (" << vals.size()
//                    << ") differs from the number of points ("
//                    <<  mapperPtr_().sourceSize()
//                    << ") in file " << vals.objectPath() << exit(FatalError);
//            }

//            startAverage_ = vals.average();

            // Reread values and interpolate
            fileName valsFile
            (
                this->db().time().constant()
               /"boundaryData"
               /this->patch().name()
               /sampleTimes_[startSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                startAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

//            // Reread values and interpolate
//            AverageIOField<Type> vals
//            (
//                IOobject
//                (
//                    fieldTableName_,
//                    this->db().time().constant(),
//                    "boundaryData"
//                   /this->patch().name()
//                   /sampleTimes_[endSampleTime_].name(),
//                    this->db(),
//                    IOobject::MUST_READ,
//                    IOobject::AUTO_WRITE,
//                    false
//                )
//            );

//            if (vals.size() != mapperPtr_().sourceSize())
//            {
//                FatalErrorInFunction
//                    << "Number of values (" << vals.size()
//                    << ") differs from the number of points ("
//                    <<  mapperPtr_().sourceSize()
//                    << ") in file " << vals.objectPath() << exit(FatalError);
//            }

//            endAverage_ = vals.average();

            // Reread values and interpolate
            fileName valsFile
            (
                this->db().time().constant()
               /"boundaryData"
               /this->patch().name()
               /sampleTimes_[endSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                endAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    scalar tend = this->db().time().elapsedCpuTime();

    Info << "myBC takes time = " << tend-tstart << endl;
}


template<class Type>
void Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    checkTable();

    // Interpolate between the sampled data

    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }


        if(setScale_)
        {
            if (meanSampledValues_.size())
                Info << "startSampledValues_last= " 
                        << startSampledValues_[meanSampledValues_.size()-1] << endl;

            startSampledValues_ = (cmptMultiply( ScaleFactorRMS_ , 
                                                 (startSampledValues_-meanSampledValues_)
                                               ) + meanSampledValues_
                                  ) * ScaleFactorMEAN_;
            if (meanSampledValues_.size())
                Info << "startSampledValues_last_after= " 
                        << startSampledValues_[meanSampledValues_.size()-1] << endl;

        }

        this->operator==(startSampledValues_);


        wantedAverage = startAverage_;
    }
    else
    {

        scalar timeVal = this->db().time().value();
        timeVal = timeVal - nRecycles_*DeltaT_;

        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (timeVal - start)/(end - start);

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        if(setScale_)
        {
            Field<Type> B4scaledValues(((1 - s)*startSampledValues_ + s*endSampledValues_));

            this->operator==( (cmptMultiply(ScaleFactorRMS_ , 
                                                   (B4scaledValues-meanSampledValues_)
                                        ) + meanSampledValues_) * ScaleFactorMEAN_
                            );
            label meanFieldsize = meanSampledValues_.size();

            if (meanFieldsize)
            {
                Pout << "meanSampledValues_size =" 
                    << meanSampledValues_.size() << endl;

                Pout << "meanSampledValues_last =" 
                    << meanSampledValues_[meanFieldsize-1] << endl;

                Pout << "B4scaled_last = "  
                    << B4scaledValues[meanFieldsize-1] << endl;

                Pout << "Afscaled_last = " 
                    << this->operator[](meanFieldsize-1) << endl;
            }

        }
        else
        {
        this->operator==((1 - s)*startSampledValues_ + s*endSampledValues_);
        }

        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const Field<Type>& fld = *this;

        Type averagePsi =
            gSum(this->patch().magSf()*fld)
           /gSum(this->patch().magSf());

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            this->operator==(fld + offset);
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            this->operator==(scale*fld);
        }
    }

//    // apply offset to mapped values
//    const scalar t = this->db().time().timeOutputValue();
//    this->operator==(*this + offset_->value(t));

//    if (debug)
//    {
//        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
//            << " max:" << gMax(*this)
//            << " avg:" << gAverage(*this) << endl;
//    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        this->operator==(*this + offset_->value(t));
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::recycledTimeVaryingMappedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);

//    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
//    if (perturb_ != 1e-5)
//    {
//        os.writeKeyword("perturb") << perturb_ << token::END_STATEMENT << nl;
//    }

//    if (fieldTableName_ != this->internalField().name())
//    {
//        os.writeKeyword("fieldTable") << fieldTableName_
//            << token::END_STATEMENT << nl;
//    }

//    if
//    (
//        (
//           !mapMethod_.empty()
//         && mapMethod_ != "planarInterpolation"
//        )
//    )
//    {
//        os.writeKeyword("mapMethod") << mapMethod_
//            << token::END_STATEMENT << nl;
//    }

//    offset_->writeData(os);

    this->writeEntryIfDifferent(os, "setAverage", Switch(false), setAverage_);

    this->writeEntryIfDifferent(os, "setScale", Switch(false), setScale_);

    if (setScale_)
    {
        os.writeKeyword("ScaleFactorRMS") << ScaleFactorRMS_ 
            << token::END_STATEMENT << nl;

        os.writeKeyword("ScaleFactorMEAN") << ScaleFactorMEAN_ 
            << token::END_STATEMENT << nl;
    }

    this->writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);

    this->writeEntryIfDifferent
    (
        os,
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );

    this->writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );

    if (offset_.valid())
    {
        offset_->writeData(os);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
