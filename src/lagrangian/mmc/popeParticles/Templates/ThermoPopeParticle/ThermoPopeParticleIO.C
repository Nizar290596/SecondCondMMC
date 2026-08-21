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

#include "ThermoPopeParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::ThermoPopeParticle<ParticleType>::propertyList_ =
    Foam::ThermoPopeParticle<ParticleType>::propertyList();

template<class ParticleType>
Foam::HashTable<Foam::label, Foam::word> Foam::ThermoPopeParticle<ParticleType>::indexInXiC_= 
    Foam::HashTable<Foam::label, Foam::word>(0);
    
template<class ParticleType>
Foam::List<Foam::word> Foam::ThermoPopeParticle<ParticleType>::XiCNames_= 
    Foam::List<Foam::word>(0);
    
template<class ParticleType>
Foam::List<Foam::word>  Foam::ThermoPopeParticle<ParticleType>::componentNames_= 
    Foam::List<Foam::word>(0);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ThermoPopeParticle<ParticleType>::ThermoPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParticleType(mesh, is, readFields, newFormat),
    T_(0.0),
    Y_(0),
    XiC_(0),
    hA_(0),
    hEqv_(0),
    dpYsource_(0),
    dpXiCsource_(0),
    dp_hAsource_(0),
    pc_(0.0),
    NumActSp_(0)
{
    if (readFields)
    {
        is >> T_
           >> Y_
           >> XiC_        
           >> hA_
           >> hEqv_
           >> dpYsource_
           >> dpXiCsource_
           >> dp_hAsource_
           >> pc_
           >> NumActSp_;
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "ThermoPopeParticle::ThermoPopeParticle(const polyMesh&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * * Read / Write * * * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::ThermoPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
}

template<class ParticleType>
template<class CloudType, class CompositionType>
void Foam::ThermoPopeParticle<ParticleType>::readFields
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
    
    IOField<scalar> T(c.newIOobject("T", IOobject::MUST_READ));
    c.checkFieldIOobject(c, T);

    IOField<scalar> hA(c.newIOobject("hA", IOobject::MUST_READ));
    c.checkFieldIOobject(c, hA);

    IOField<scalar> hEqv(c.newIOobject("hEqv", IOobject::MUST_READ));
    c.checkFieldIOobject(c, hEqv);

//    label numfV = c.numfVsoot();

//    List<IOField<scalar>> fVsoot(numfV);
//    
//    forAll(fVsoot, ii)
//    {
//        fVsoot[ii] = c.newIOobject(word("fVsoot"+ii), IOobject::MUST_READ);
//        c.checkFieldIOobject(c, fVsoot[ii]);
//    }
    
    IOField<scalar> pc(c.newIOobject("pc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, pc);

    IOField<scalar> NumActSp(c.newIOobject("NumActSp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, NumActSp);

    label i = 0;
    forAllIters(c, iter)
    {
        ThermoPopeParticle<ParticleType>& p = iter();

        p.T_            = T[i];
        p.hA_           = hA[i];
        p.hEqv_         = hEqv[i];
//        p.fVsoot_       = fVsoot[i];
        p.dp_hAsource_  = 0.0;
        p.pc_           = pc[i];
        p.NumActSp_     = NumActSp[i];
        
        i++;
    }


    // Get names and sizes for each Y...
    const wordList& compTypes = compModel.componentNames();

    const label nComp = compTypes.size();

    // Set storage for each Y... for each particle
    forAllIters(c, iter)
    {
        ThermoPopeParticle<ParticleType>& p = iter();
        p.Y_.setSize(nComp, 0.0);
        p.dpYsource_.setSize(nComp, 0.0);
    }

    // Populate Y for each parcel
    forAll(compTypes, j)
    {
        IOField<scalar> Y
        (
            c.newIOobject
            (
                "Y" + compTypes[j],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            ThermoPopeParticle<ParticleType>& p = iter();
            
            p.Y_[j] = Y[i];
            
            p.dpYsource_[j] = 0.0;//dpYsource[i];
            
            i++;
        }
    }
    
    // Get names and sizes for each XiC...
    const wordList& XiCTypes = c.coupling().XiCNames();

    const label noXiC = XiCTypes.size();

    // Set storage for each XiC... for each particle
    forAllIters(c, iter)
    {
        ThermoPopeParticle<ParticleType>& p = iter();
        p.XiC_.setSize(noXiC, 0.0);
        p.dpXiCsource_.setSize(noXiC, 0.0);
    }

    // Populate XiC and dXiC for each parcel !!!

    forAllConstIters(XiCTypes, XiC_i)
    {
        IOField<scalar> XiC
        (
            c.newIOobject
            (
                *XiC_i,
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            ThermoPopeParticle<ParticleType>& p = iter();
            
            p.XiC(*XiC_i) = XiC[i];
            p.dpXiCsource(*XiC_i) = 0.0;
            
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
void Foam::ThermoPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


template<class ParticleType>
template<class CloudType, class CompositionType>
void Foam::ThermoPopeParticle<ParticleType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    //ParticleType::writeFields(c);

    label np = c.size();
    if (np > 0)
    {
        IOField<scalar> T(c.newIOobject("T", IOobject::NO_READ), np);
        IOField<scalar> hA(c.newIOobject("hA", IOobject::NO_READ), np);
        IOField<scalar> hEqv(c.newIOobject("hEqv", IOobject::NO_READ), np);
        IOField<scalar> pc(c.newIOobject("pc", IOobject::NO_READ), np);
        IOField<scalar> rho(c.newIOobject("rho", IOobject::NO_READ), np);
        IOField<scalar> NumActSp(c.newIOobject("NumActSp", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIters( c, iter)
        {
            const ThermoPopeParticle<ParticleType>& p = iter();

            T[i]      = p.T_;
            hA[i]     = p.hA_;
            hEqv[i]   = p.hEqv_;
            pc[i]     = p.pc_;
            rho[i]    = compModel.particleMixture(p.Y_).rho(p.pc_,p.T_);
            NumActSp[i] = p.NumActSp_;

            i++;
        }

        T.write();
        hA.write();
        hEqv.write();
        pc.write();
        rho.write();
        NumActSp.write();

        // Write the composition mass fractions
        const wordList& compTypes = compModel.componentNames();

        forAll(compTypes, j)
        {
            IOField<scalar> Y
            (
                c.newIOobject
                (
                    "Y" + compTypes[j],
                    IOobject::NO_READ
                ),
                np
            );
      
            label i = 0;
            forAllConstIters
            (
                c,
                iter
            )
            {
                const ThermoPopeParticle<ParticleType>& p = iter();                
                Y[i] = p.Y()[j];
                i++;
            }
            
            Y.write();

        }
        
        // Write the coupling variables and sources 
        const wordList& XiCTypes = c.coupling().XiCNames();

        forAllConstIter(wordList,XiCTypes, XiC_i)
        {
            IOField<scalar> XiC
            (
                c.newIOobject
                (
                    *XiC_i,
                    IOobject::NO_READ
                ),
                np
            );
      
            // If only the thermo cloud is constructed the particles 
            // have not called setStaticProperties
            if (indexInXiC_.size() == 0)
                setStaticProperties(c);
      
            label i = 0;
            forAllConstIters
            (
                c,
                iter
            )
            {
                const ThermoPopeParticle<ParticleType>& p = iter();            
                XiC[i] = p.XiC(*XiC_i);
                i++;
            }
            
            XiC.write();
            
        }
        
        
        // Write the composition mole fractions if required
        if (!compModel.printMoleFractionsEnabled())
            return;

        forAll(compTypes, j)
        {
            IOField<scalar> X
            (
                c.newIOobject
                (
                    "X" + compTypes[j],
                    IOobject::NO_READ
                ),
                np
            );
      
            label i = 0;
            forAllConstIters
            (
                c,
                iter
            )
            {
                const ThermoPopeParticle<ParticleType>& p = iter();                
                X[i] = compModel.X(p.Y())[j];
                i++;
            }
            
            X.write();

        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ThermoPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.T()
        << token::SPACE << p.Y()
        << token::SPACE << p.XiC()        
        << token::SPACE << p.hA()
        << token::SPACE << p.hEqv()
        << token::SPACE << p.dpYsource()
        << token::SPACE << p.dpXiCsource()
        << token::SPACE << p.dp_hAsource()
        << token::SPACE << p.pc()
        << token::SPACE << p.NumActSp();

    // Check state of Ostream
    os.check
    (
       "Ostream& operator<<(Ostream&, const ThermoPopeParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //
