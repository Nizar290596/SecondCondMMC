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

#include "MixingPopeParticle.H"
#include "IOstreams.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::MixingPopeParticle<ParticleType>::propertyList_ =
    Foam::MixingPopeParticle<ParticleType>::propertyList();


template<class ParticleType>
Foam::HashTable<Foam::label, Foam::word> Foam::MixingPopeParticle<ParticleType>::indexInXiR_= 
    Foam::HashTable<Foam::label, Foam::word>(0);
    
template<class ParticleType>
Foam::List<Foam::word> Foam::MixingPopeParticle<ParticleType>::XiRNames_=
    Foam::List<Foam::word>(0);

//template<class ParticleType>
//Foam::scalar Foam::MixingPopeParticle<ParticleType>::secondCondBeta_s_ = 1.0;
    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::MixingPopeParticle<ParticleType>::MixingPopeParticle
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
    secondCondFlag_(0),
    omegaOU_(0.0),
    phi_(0.0),
    phiModified_(0.0),
    burnedAge_(0.0)
{
    if (readFields)
    {
        is  >>  XiR_
            >> dXiR_
            >> dx_
	    >> secondCondFlag_
            >> omegaOU_
	    >> phi_
            >> phiModified_
            >> burnedAge_;
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "MixingPopeParticle::MixingPopeParticle(const polyMesh&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * * Read / Write * * * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::MixingPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
}


template<class ParticleType>
template<class CloudType, class MixingType>
void Foam::MixingPopeParticle<ParticleType>::readFields
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

    IOField<scalar> dx(c.newIOobject("dx", IOobject::READ_IF_PRESENT), scalarField(c.size(), 0.0));
    c.checkFieldIOobject(c, dx);
    const dictionary sc = c.cloudProperties().subOrEmptyDict("secondConditioning");
    const word restartMode = sc.lookupOrDefault<word>("restartMode", "resume");
    if (restartMode != "resume" && restartMode != "initialize")
        FatalErrorInFunction << "secondConditioning/restartMode must be resume or initialize" << exit(FatalError);
    const bool initialize = c.secondCondMixingEnabled() && restartMode == "initialize";
    const bool mustRead = c.secondCondMixingEnabled() && !initialize;
    IOField<label> secondCondFlag(c.newIOobject("secondCondFlag",
        mustRead ? IOobject::MUST_READ : IOobject::READ_IF_PRESENT), Field<label>(c.size(), label(0)));
    IOField<scalar> omegaOU(c.newIOobject("omegaOU",
        mustRead ? IOobject::MUST_READ : IOobject::READ_IF_PRESENT), scalarField(c.size(), 0.0));
    IOField<scalar> phi(c.newIOobject("phi",
        mustRead ? IOobject::MUST_READ : IOobject::READ_IF_PRESENT), scalarField(c.size(), 0.0));
    IOField<scalar> phiModified(c.newIOobject("phiModified",
        mustRead ? IOobject::MUST_READ : IOobject::READ_IF_PRESENT), scalarField(c.size(), 0.0));
    IOField<scalar> burnedAge(c.newIOobject("burnedAge",
        mustRead && sc.lookupOrDefault<scalar>("burnedAgeRate", 0) > 0
            ? IOobject::MUST_READ : IOobject::READ_IF_PRESENT), scalarField(c.size(), 0.0));
    c.checkFieldIOobject(c, secondCondFlag);
    c.checkFieldIOobject(c, omegaOU);
    c.checkFieldIOobject(c, phi);
    c.checkFieldIOobject(c, phiModified);
    c.checkFieldIOobject(c, burnedAge);
    label i = 0;
    forAllIters(c, iter)
    {
        MixingPopeParticle<ParticleType>& p = iter();
        p.dx_ = dx[i];
        p.secondCondFlag_ = secondCondFlag[i];
        p.omegaOU_ = omegaOU[i];
        p.phi_ = phi[i];
        p.phiModified_ = phiModified[i];
        p.burnedAge_ = burnedAge[i];
        if (initialize) c.initializeSecondConditioningState(iter());
        if (c.secondCondMixingEnabled() &&
            (!(p.secondCondFlag_ == 0 || p.secondCondFlag_ == 1)
             || !(p.phi_ >= 0 && p.phi_ <= 1 && p.phiModified_ >= 0 && p.burnedAge_ >= 0)
             || !std::isfinite(p.omegaOU_ + p.phiModified_ + p.burnedAge_)))
            FatalErrorInFunction << "Invalid second-conditioning restart state" << exit(FatalError);
        ++i;
    }
    if (initialize)
        Info<< "Initialized second-conditioning reference state and subset from restart temperature. "
            << "Use restartMode resume for subsequent restarts." << nl;

    // Get names and sizes for each XiR...
    const wordList& XiRTypes = mixModel.XiRNames();

    const label noXiR = XiRTypes.size();
    
    const HashTable<label, word>& XiRIndexes = mixModel.XiR().rVarInXiR();

    // Set storage for each XiR... for each particle
    forAllIters(c, iter)
    {
        MixingPopeParticle<ParticleType>& p = iter();
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
            MixingPopeParticle<ParticleType>& p = iter();
            
            p.XiR_[XiRIndexes[*XiR_i]] = XiR[i];
            p.dXiR_[XiRIndexes[*XiR_i]] = dXiR[i];
            
            i++;
        }
    }

    // Initialize the particle sampling
    for (auto& p : c)
    {
        p.initStatisticalSampling();
    }
}
    

template<class ParticleType>
template<class CloudType>
void Foam::MixingPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


template<class ParticleType>
template<class CloudType, class MixingType>
void Foam::MixingPopeParticle<ParticleType>::writeFields
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
	IOField<label>  secondCondFlag
        (
            c.newIOobject("secondCondFlag", IOobject::NO_READ), np
        );
        IOField<scalar> omegaOU(c.newIOobject("omegaOU", IOobject::NO_READ), np);
	IOField<scalar> phi(c.newIOobject("phi", IOobject::NO_READ), np);
        IOField<scalar> phiModified
        (
            c.newIOobject("phiModified", IOobject::NO_READ), np
        );

        IOField<scalar> burnedAge(c.newIOobject("burnedAge", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIters(c, iter)
        {
            const MixingPopeParticle<ParticleType>& p = iter();
            dx[i]             = p.dx_;
	    secondCondFlag[i] = p.secondCondFlag_;
            omegaOU[i]        = p.omegaOU_;
	    phi[i]            = p.phi_;
            phiModified[i]    = p.phiModified_;
            burnedAge[i]       = p.burnedAge_;
            i++;
        }

        dx.write();
	secondCondFlag.write();
        omegaOU.write();
	phi.write();
        phiModified.write();
        burnedAge.write();
        
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
                const MixingPopeParticle<ParticleType>& p = iter();
                
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
    const MixingPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.XiR_
        << token::SPACE << p.dXiR_
        << token::SPACE << p.dx_
	<< token::SPACE << p.secondCondFlag_
        << token::SPACE << p.omegaOU_
	<< token::SPACE << p.phi_
        << token::SPACE << p.phiModified_
        << token::SPACE << p.burnedAge_;
 
    // Check state of Ostream
    os.check
    (
       "Ostream& operator<<(Ostream&, const MixingPopeParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //

