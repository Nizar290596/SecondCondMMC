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

#include "SootPopeParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::SootPopeParticle<ParticleType>::SootPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParticleType(mesh, is, readFields, newFormat),
    fVsoot_(0),
    Intermittency_(0),
    sootVf_(0.),
    sootIm_(1.),
    ySoot_(0.),
    nSoot_(0.),
    sootContribs_(6, 0.),
    XO2XN2_(0.)
{
    if (readFields)
    {
        is >> fVsoot_
           >> Intermittency_
           >> sootVf_
           >> sootIm_
           >> ySoot_
           >> nSoot_
           >> sootContribs_
           >> XO2XN2_;
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "SootPopeParticle::SootPopeParticle(const polyMesh&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * * Read / Write * * * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::SootPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
}

template<class ParticleType>
template<class CloudType, class CompositionType>
void Foam::SootPopeParticle<ParticleType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
    
    setStaticProperties(c);

    IOField<scalar> sootVf(c.newIOobject("sootVf", IOobject::MUST_READ));
    c.checkFieldIOobject(c, sootVf);

    IOField<scalar> sootIm(c.newIOobject("sootIm", IOobject::MUST_READ));
    c.checkFieldIOobject(c, sootIm);

    IOField<scalar> ySoot(c.newIOobject("ySoot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ySoot);

    IOField<scalar> nSoot(c.newIOobject("nSoot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, nSoot);

    IOField<scalar> XO2XN2(c.newIOobject("XO2XN2", IOobject::MUST_READ));
    c.checkFieldIOobject(c, XO2XN2);

    label i = 0;
    forAllIters(c, iter)
    {
        auto& p = iter();

        p.sootVf_       = sootVf[i];
        p.sootIm_       = sootIm[i];
        p.ySoot_        = ySoot[i];
        p.nSoot_        = nSoot[i];
        p.XO2XN2_       = XO2XN2[i];

        i++;
    }



    //Get sizes of fVsoot_
    label numfV = c.numfVsoot();

    forAllIters(c, iter)
    {
        SootPopeParticle<ParticleType>& p = iter();
        p.fVsoot_.setSize(numfV, 0.0);
    }

    for(label ii=0; ii<numfV; ii++)
    {
        IOField<scalar> fVsoot
        (
            c.newIOobject
            (
                "fVsoot" + std::to_string(ii),
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            SootPopeParticle<ParticleType>& p = iter();
            
            p.fVsoot_[ii] = fVsoot[i];
            
            i++;
        }

    }

    // Get sizes of Intermittency_

    forAllIters(c, iter)
    {
        SootPopeParticle<ParticleType>& p = iter();
        p.Intermittency_.setSize(numfV, 0.0);
    }

    for(label ii=0; ii<numfV; ii++)
    {
        IOField<scalar> Intermittency
        (
            c.newIOobject
            (
                "Intermittency" + std::to_string(ii),
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            SootPopeParticle<ParticleType>& p = iter();
            
            p.Intermittency_[ii] = Intermittency[i];
            
            i++;
        }

    }

    // Initialize the particle sampling
    for (auto p : c)
    {
        p.initStatisticalSampling();
    }
}


template<class ParticleType>
template<class CloudType>
void Foam::SootPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


template<class ParticleType>
template<class CloudType, class CompositionType>
void Foam::SootPopeParticle<ParticleType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParticleType::writeFields(c,compModel);

    label np = c.size();
    if (np > 0)
    {
        IOField<scalar> sootVf(c.newIOobject("sootVf", IOobject::NO_READ), np);
        IOField<scalar> sootIm(c.newIOobject("sootIm", IOobject::NO_READ), np);
        IOField<scalar> ySoot(c.newIOobject("ySoot", IOobject::NO_READ), np);
        IOField<scalar> nSoot(c.newIOobject("nSoot", IOobject::NO_READ), np);
        IOField<scalar> XO2XN2(c.newIOobject("XO2XN2", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIters(c, iter)
        {
            const auto& p = iter();

            sootVf[i] = p.sootVf_;
            sootIm[i] = p.sootIm_;
            ySoot[i]  = p.ySoot_;
            nSoot[i]  = p.nSoot_;
            XO2XN2[i] = p.XO2XN2_;

            i++;
        }

        sootVf.write();
        sootIm.write();
        ySoot.write();
        nSoot.write();
        XO2XN2.write();

        // write the soot contributions from each step
        List<string> sootContribsTag = 
            {"yInC2H2", "ySgC2H2", "yOxO2", "yOxOH", "nInC2H2", "nCg"};

        for(int j = 0; j < sootContribsTag.size(); j++)
        {
            IOField<scalar> sootContribs
            (
                c.newIOobject
                (
                    "sootContribs_" + sootContribsTag[j],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<SootPopeParticle<ParticleType> >,
                c,
                iter
            )
            {
                const SootPopeParticle<ParticleType>& p = iter();
                sootContribs[i] = p.sootContribs_[j];
                i++;
            }
            sootContribs.write();
        }

        // Write soot volume fraction
        const label numfV = c.numfVsoot();

        for(label ii=0; ii<numfV; ii++)
        {
            IOField<scalar> fVsoot
            (
                c.newIOobject
                (
                    "fVsoot" + std::to_string(ii),
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIters(c, iter)
            {
                const SootPopeParticle<ParticleType>& p = iter();                
                fVsoot[i] = p.fVsoot_[ii];
                i++;
            }
            
            fVsoot.write();
        }

        //write Intermittency
        for(label ii=0; ii<numfV; ii++)
        {
            IOField<scalar> Intermittency
            (
                c.newIOobject
                (
                    "Intermittency" + std::to_string(ii),
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIters(c, iter)
            {
                const SootPopeParticle<ParticleType>& p = iter();                
                Intermittency[i] = p.Intermittency_[ii];
                i++;
            }
            
            Intermittency.write();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SootPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.fVsoot()
        << token::SPACE << p.Intermittency()
        << token::SPACE << p.sootVf()
        << token::SPACE << p.sootIm()
        << token::SPACE << p.ySoot()
        << token::SPACE << p.nSoot()
        << token::SPACE << p.sootContribs()
        << token::SPACE << p.XO2XN2();

    // Check state of Ostream
    os.check
    (
       "Ostream& operator<<(Ostream&, const SootPopeParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //
