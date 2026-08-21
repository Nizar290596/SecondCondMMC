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

#include "PremixedMixingPopeParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::PremixedMixingPopeParticle<ParticleType>::propertyList_ =
    Foam::PremixedMixingPopeParticle<ParticleType>::propertyList();


template<class ParticleType>
Foam::HashTable<Foam::label, Foam::word> Foam::PremixedMixingPopeParticle<ParticleType>::indexInXiR_= 
    Foam::HashTable<Foam::label, Foam::word>(0);
    
template<class ParticleType>
Foam::List<Foam::word> Foam::PremixedMixingPopeParticle<ParticleType>::XiRNames_= 
    Foam::List<Foam::word>(0);
    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::PremixedMixingPopeParticle<ParticleType>::PremixedMixingPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParticleType(mesh, is, readFields, newFormat),
    XiR_(0),
    dXiR_(0),
    dx_(0.0),
    mixTime_(0.0),
    mixExt_(0.0),
    progVarInt_(0.0),
    XiRrlxFactor_(0.0),
    rlxStep_(0)
{
    if (readFields)
    {
        is  >>  XiR_
            >> dXiR_
            >> dx_
            >> mixTime_
            >> mixExt_
            >> progVarInt_
            >> XiRrlxFactor_
            >> rlxStep_;
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "PremixedMixingPopeParticle::PremixedMixingPopeParticle(const polyMesh&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * * Read / Write * * * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::PremixedMixingPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
}


template<class ParticleType>
template<class CloudType, class MixingType>
void Foam::PremixedMixingPopeParticle<ParticleType>::readFields
(
    CloudType& c,
    const MixingType& mixModel
)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
    
    setStaticProperties(c);

    IOField<scalar> dx(c.newIOobject("dx", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dx);

    IOField<scalar> mixTime(c.newIOobject("mixTime", IOobject::MUST_READ));
    c.checkFieldIOobject(c, mixTime);

    IOField<scalar> mixExt(c.newIOobject("mixExt", IOobject::MUST_READ));
    c.checkFieldIOobject(c, mixExt);

    IOField<scalar> progVarInt(c.newIOobject("progVarInt", IOobject::MUST_READ));
    c.checkFieldIOobject(c, progVarInt);

    IOField<scalar> XiRrlxFactor(c.newIOobject("XiRrlxFactor", IOobject::MUST_READ));
    c.checkFieldIOobject(c, XiRrlxFactor);

    IOField<scalar> rlxStep(c.newIOobject("rlxStep", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rlxStep);

    label i = 0;
    forAllIters(c, iter)
    {
        PremixedMixingPopeParticle<ParticleType>& p = iter();
        p.dx_           = dx[i];
        p.mixTime_      = mixTime[i];
        p.mixExt_       = mixExt[i];
        p.progVarInt_   = progVarInt[i];
        p.XiRrlxFactor_ = XiRrlxFactor[i];
        p.rlxStep_      = rlxStep[i];       
        i++;
    }
    
    // Get names and sizes for each XiR...
    const wordList& XiRTypes = mixModel.XiRNames();

    const label noXiR = XiRTypes.size();
    
    const HashTable<label, word>& XiRIndexes = mixModel.XiR().rVarInXiR();

    // Set storage for each XiR... for each particle
    forAllIters(c, iter)
    {
        PremixedMixingPopeParticle<ParticleType>& p = iter();
        p.XiR_.setSize(noXiR, 0.0);
        p.dXiR_.setSize(noXiR, 0.0);
    }

    // Populate XiR and dXiR for each parcel !!!
    forAllConstIter(wordList,XiRTypes, XiR_i)
    {
        IOField<scalar> XiR
        (
            c.newIOobject
            (
                *XiR_i,
                IOobject::MUST_READ
            )
        );

        IOField<scalar> dXiR
        (
            c.newIOobject
            (
                "d" + *XiR_i,
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            PremixedMixingPopeParticle<ParticleType>& p = iter();
            
            p.XiR_[XiRIndexes[*XiR_i]] = XiR[i];
            p.dXiR_[XiRIndexes[*XiR_i]] = dXiR[i];
            
            i++;
        }
    }
    
}
    

template<class ParticleType>
template<class CloudType>
void Foam::PremixedMixingPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


template<class ParticleType>
template<class CloudType, class MixingType>
void Foam::PremixedMixingPopeParticle<ParticleType>::writeFields
(
    const CloudType& c,
    const MixingType& mixModel
)
{
    ParticleType::writeFields(c);

    label np = c.size();
    if (np > 0)
    {
        IOField<scalar> dx(c.newIOobject("dx", IOobject::NO_READ), np);
        IOField<scalar> mixTime(c.newIOobject("mixTime", IOobject::NO_READ), np);
        IOField<scalar> mixExt(c.newIOobject("mixExt", IOobject::NO_READ), np);
        IOField<scalar> progVarInt(c.newIOobject("progVarInt", IOobject::NO_READ), np);
        IOField<scalar> XiRrlxFactor(c.newIOobject("XiRrlxFactor", IOobject::NO_READ), np);
        IOField<scalar> rlxStep(c.newIOobject("rlxStep", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIters(c, iter)
        {
            const PremixedMixingPopeParticle<ParticleType>& p = iter();
            dx[i] = p.dx_;
            mixTime[i] = p.mixTime_;
            mixExt[i] = p.mixExt_;
            progVarInt[i] = p.progVarInt_;
            XiRrlxFactor[i] = p.XiRrlxFactor_;
            rlxStep[i] = p.rlxStep_;
            i++;
        }

        dx.write();
        mixTime.write();
        mixExt.write();
        progVarInt.write();
        XiRrlxFactor.write();
        rlxStep.write();
        
        // Write the reference variables and distances in Xi space
        const wordList& XiRTypes = mixModel.XiRNames();
        
        const HashTable<label, word>& XiRIndexes = mixModel.XiR().rVarInXiR();

        forAllConstIter(wordList,XiRTypes, XiR_i)
        {
            IOField<scalar> XiR
            (
                c.newIOobject
                (
                    *XiR_i,
                    IOobject::NO_READ
                ),
                np
            );

            IOField<scalar> dXiR
            (
                c.newIOobject
                (
                    "d" + *XiR_i,
                    IOobject::NO_READ
                ),
                np
            );
      
            label k = 0;
            forAllConstIters
            (
                c,
                iter
            )
            {
                const PremixedMixingPopeParticle<ParticleType>& p = iter();
                
                XiR[k]  = p.XiR()[XiRIndexes[*XiR_i]];                
                dXiR[k] = p.dXiR()[XiRIndexes[*XiR_i]];

                k++;
            }
            
            XiR.write();
            
            dXiR.write();
        }
        
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PremixedMixingPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.XiR_
        << token::SPACE << p.dXiR_
        << token::SPACE << p.dx_
        << token::SPACE << p.mixTime_
        << token::SPACE << p.mixExt_
        << token::SPACE << p.progVarInt_
        << token::SPACE << p.XiRrlxFactor_
        << token::SPACE << p.rlxStep_;
 
    // Check state of Ostream
    os.check
    (
       "Ostream& operator<<(Ostream&, const PremixedMixingPopeParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //
