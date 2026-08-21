/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    popeParticle

Description

 ICE Revision: $Id: popeCloudIO.C 7092 2007-01-25 21:38:27Z bgschaid $ 
\*---------------------------------------------------------------------------*/

#include "ItoPopeParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from Istream
template<class ParticleType>
Foam::ItoPopeParticle<ParticleType>::ItoPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParticleType(mesh, is,readFields,newFormat),

    sCell_(0),

    mass_(0.0),
    
    dpMsource_(0.0),

    wt_(0.0),
    
    A_(vector::zero),

    rw_(vector::zero),

    toPos_(vector::zero)
{
    if (readFields)
    {
        is  >> sCell_
            >> mass_
            >> dpMsource_
            >> wt_
            >> A_
            >> rw_
            >> toPos_;
    }

    initStatisticalSampling();
    
    is.check("ItoPopeParticle::ItoPopeParticle(IStream &");
}


template<class ParticleType>
template<class CloudType> 
void Foam::ItoPopeParticle<ParticleType>::writeFields(const CloudType &c)
{
    ParticleType::writeFields(c);

    label np = c.size();

    IOField<scalar> sCell(c.newIOobject("sCell",IOobject::NO_READ),np);
    IOField<scalar> m(c.newIOobject("m",IOobject::NO_READ),np);
    IOField<scalar> wt(c.newIOobject("wt",IOobject::NO_READ),np);
    IOField<vector> A(c.newIOobject("A",IOobject::NO_READ),np);
    IOField<vector> rw(c.newIOobject("rw",IOobject::NO_READ),np);
    IOField<vector> toPos(c.newIOobject("toPos",IOobject::NO_READ),np);

    label i = 0;

    forAllConstIters(c,iter)
    {
        const ItoPopeParticle<ParticleType>& p = iter();

        sCell[i] = p.sCell();

        m[i] = p.m();

        wt[i] = p.wt();
        
        A[i] = p.A();
             
        rw[i] = p.rw();
        
        toPos[i] = p.toPos();

        i++;
    }

    sCell.write();

    m.write();

    wt.write();
    
    A.write();

    rw.write();

    toPos.write();
}


template<class ParticleType>
template<class CloudType>
void Foam::ItoPopeParticle<ParticleType>::readFields(CloudType &c)
{
    if(!c.size())
    {
        return;
    }
    
    ParticleType::readFields(c);

    IOField<scalar> sCell(c.newIOobject("sCell",IOobject::MUST_READ));
    c.checkFieldIOobject(c, sCell);

    IOField<scalar> m(c.newIOobject("m",IOobject::MUST_READ));
    c.checkFieldIOobject(c, m);

    IOField<scalar> wt(c.newIOobject("wt",IOobject::MUST_READ));
    c.checkFieldIOobject(c, wt);

    IOField<vector> A(c.newIOobject("A",IOobject::MUST_READ));
    c.checkFieldIOobject(c, A);

    IOField<vector> rw(c.newIOobject("rw",IOobject::MUST_READ));
    c.checkFieldIOobject(c, rw);

    IOField<vector> toPos(c.newIOobject("toPos",IOobject::MUST_READ));
    c.checkFieldIOobject(c, toPos);

    label i = 0;
    forAllIters(c,iter)
    {
        ItoPopeParticle<ParticleType>& p = iter();

        p.sCell() = sCell[i];

        p.m() = m[i];
        
        p.dpMsource() = 0.0;//dpMsource[i];

        p.wt() = wt[i];
        
        p.A() = A[i];

        p.rw() = rw[i];

        p.toPos() = toPos[i];

        i++;
    }

    // Initialize the particle sampling
    for (auto p : c)
    {
        p.initStatisticalSampling();
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os, 
    const ItoPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.sCell()
        << token::SPACE << p.m()
        << token::SPACE << p.dpMsource()
        << token::SPACE << p.wt()
        << token::SPACE << p.A()
        << token::SPACE << p.rw()
        << token::SPACE << p.toPos();

    // Check state of Ostream
    os.check("operator<<(Ostream& os, const popeParticle& p)");
        
    return os;
}


// ************************************************************************* //
