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

#include "sprayThermoParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::sprayThermoParcel<ParcelType>::propertyList_ =
    Foam::sprayThermoParcel<ParcelType>::propertyList();

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields),
    TD_(0.0),
    TGas_(0.0),
    TFilm_(0.0),
    TSurf_(0.0),
    pc_(0.0),
    fGas_(0.0),
    fSurf_(0.0),
    fFilm_(0.0),
    nThermoSpecies_(0),
    nFuelSpecies_(0),
    YFuel_(3),
    YFuelVap_(3),
    YGas_(3),
    YFilm_(3),
    YSurf_(3),
    haGas_(0.0),
    hD_(0.0),
    haSurf_(0.0),
    haFilm_(0.0),
    haFuelVap_(0.0),
    muFilm_(0.0),
    rhoFilm_(0.0),
    qR_(0.0),
    qD_(0.0),
    mFlux_(0.0),
    molWtSstate_(0.0),
    ZFilm_(0.0),
    BM_(0.0) 
{
    if (readFields)
    {

        TD_ = readScalar(is);
        TGas_ = readScalar(is);
        TFilm_ = readScalar(is);
        TSurf_ = readScalar(is);
        pc_ = readScalar(is);
        fGas_ = readScalar(is);
        fSurf_ = readScalar(is);
        fFilm_ = readScalar(is);
        nThermoSpecies_ = readLabel(is);
	    nFuelSpecies_ = readLabel(is);
        YFuel_.setSize(nFuelSpecies_,0.0);
            is >> YFuel_;
        haGas_ = readScalar(is);
        hD_ = readScalar(is);
        haSurf_ = readScalar(is);
        haFilm_ = readScalar(is);
        qR_ = readScalar(is);
        qD_ = readScalar(is);
        mFlux_ = readScalar(is);
        molWtSstate_ = readScalar(is);
        BM_ = readScalar(is);
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "sprayThermoParcel::sprayThermoParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::sprayThermoParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::sprayThermoParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    if (!c.size())
    {
        return;
    }

    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> TD(c.newIOobject("TD", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, TD);

    IOField<scalar> TGas(c.newIOobject("TGas", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, TGas);

    IOField<scalar> TFilm(c.newIOobject("TFilm", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, TFilm);

    IOField<scalar> TSurf(c.newIOobject("TSurf", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, TSurf);

    IOField<scalar> pc(c.newIOobject("pc", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, pc);

    IOField<scalar> fGas(c.newIOobject("fGas", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, fGas);

    IOField<scalar> fSurf(c.newIOobject("fSurf", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, fSurf);

    IOField<scalar> fFilm(c.newIOobject("fFilm", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, fFilm);

    IOField<label> nThermoSpecies(c.newIOobject("nThermoSpecies", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, nThermoSpecies);
    
    IOField<label> nFuelSpecies(c.newIOobject("nFuelSpecies", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, nFuelSpecies);

    IOField<scalar> haGas(c.newIOobject("haGas", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, haGas);

    IOField<scalar> hD(c.newIOobject("hD", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, hD);

    IOField<scalar> haSurf(c.newIOobject("haSurf", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, haSurf);

    IOField<scalar> haFilm(c.newIOobject("haFilm", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, haFilm);

    IOField<scalar> qR(c.newIOobject("qR", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, qR);

    IOField<scalar> qD(c.newIOobject("qD", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, qD);

    IOField<scalar> mFlux(c.newIOobject("mFlux", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mFlux);

    IOField<scalar> molWtSstate(c.newIOobject("molWtSstate", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, molWtSstate);

    IOField<scalar> BM(c.newIOobject("BM", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, BM);

    label i = 0;
    forAllIters( c, iter)
    {
        sprayThermoParcel<ParcelType>& p = iter();

        p.TD_ = TD[i];
        p.TGas_ = TGas[i];
        p.TFilm_ = TFilm[i];
        p.TSurf_ = TSurf[i];
        p.pc_ = pc[i];
        p.fGas_ = fGas[i];
        p.fSurf_ = fSurf[i];
        p.fFilm_ = fFilm[i];
        p.nThermoSpecies_ = nThermoSpecies[i];
        p.nFuelSpecies_ = nFuelSpecies[i];
	    p.haGas_ = haGas[i];
        p.hD_ = hD[i];
        p.haSurf_ = haSurf[i];
        p.haFilm_ = haFilm[i];
        p.qR_ = qR[i];
        p.qD_ = qD[i];
        p.mFlux_ = mFlux[i];
        p.molWtSstate_= molWtSstate[i];
        p.BM_ = BM[i]; 

        i++;
    }

    // Get names and sizes for each Y...
    const wordList& compTypes = compModel.componentNames();

    const label nComp = compTypes.size();

    // Set storage for each Y... for each Parcel
    forAllIters(c, iter)
    {
        sprayThermoParcel<ParcelType>& p = iter();
        p.YFuel_.setSize(c.nFuelSpecies(), 0.0);
        p.YFuelVap_.setSize(c.nFuelSpecies(), 0.0);
        p.YGas_.setSize(nComp, 0.0);
        p.YSurf_.setSize(nComp, 0.0);
        p.YFilm_.setSize(nComp, 0.0);
    }

    forAll(c.linkFG(), j)
    {
        IOField<scalar> YFuel
        (
            c.newIOobject
            (
                "YFuel" +  compTypes[c.linkFG()[j]],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            sprayThermoParcel<ParcelType>& p = iter();
            p.YFuel_[j] = YFuel[i];
            i++;
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::sprayThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::sprayThermoParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<scalar> TD(c.newIOobject("TD", IOobject::NO_READ), np);
    IOField<scalar> TGas(c.newIOobject("TGas", IOobject::NO_READ), np);   
    IOField<scalar> TFilm(c.newIOobject("TFilm", IOobject::NO_READ), np);   
    IOField<scalar> TSurf(c.newIOobject("TSurf", IOobject::NO_READ), np);
    IOField<scalar> pc(c.newIOobject("pc", IOobject::NO_READ), np);    
    IOField<scalar> fGas(c.newIOobject("fGas", IOobject::NO_READ), np);    
    IOField<scalar> fSurf(c.newIOobject("fSurf", IOobject::NO_READ), np);    
    IOField<scalar> fFilm(c.newIOobject("fFilm", IOobject::NO_READ), np);
    IOField<label> nThermoSpecies(c.newIOobject("nThermoSpecies", IOobject::NO_READ), np);   
    IOField<label> nFuelSpecies(c.newIOobject("nFuelSpecies", IOobject::NO_READ), np); 
    IOField<scalar> haGas(c.newIOobject("haGas", IOobject::NO_READ), np);    
    IOField<scalar> hD(c.newIOobject("hD", IOobject::NO_READ), np);   
    IOField<scalar> haSurf(c.newIOobject("haSurf", IOobject::NO_READ), np);    
    IOField<scalar> haFilm(c.newIOobject("haFilm", IOobject::NO_READ), np);   
    IOField<scalar> qR(c.newIOobject("qR", IOobject::NO_READ), np);     
    IOField<scalar> qD(c.newIOobject("qD", IOobject::NO_READ), np);  
    IOField<scalar> mFlux(c.newIOobject("mFlux", IOobject::NO_READ), np); 
    IOField<scalar> molWtSstate(c.newIOobject("molWtSstate", IOobject::NO_READ), np);
    IOField<scalar> BM(c.newIOobject("BM", IOobject::NO_READ), np); 

    label i = 0;
    forAllConstIters(c, iter)
    {
        const sprayThermoParcel<ParcelType>& p = iter();

        TD[i] = p.TD_; 
        TGas[i] = p.TGas_;
        TFilm[i] = p.TFilm_;
        TSurf[i] = p.TSurf_;
        pc[i] = p.pc_;
        fGas[i] = p.fGas_;
        fSurf[i] = p.fSurf_;
        fFilm[i] = p.fFilm_;
        nThermoSpecies[i] = p.nThermoSpecies_;
	    nFuelSpecies[i] = p.nFuelSpecies_;
        haGas[i] = p.haGas_;
        hD[i] = p.hD_;
        haSurf[i] = p.haSurf_;
        haFilm[i] = p.haFilm_;
        qR[i] = p.qR_;
        qD[i] = p.qD_;
        mFlux[i] = p.mFlux_;
        molWtSstate[i] =p.molWtSstate_;
        BM[i] = p.BM_;  
        i++;
    }

    TD.write(np > 0);
    TGas.write(np > 0);
    TFilm.write(np > 0);
    TSurf.write(np > 0);
    pc.write(np > 0);
    fGas.write(np > 0);
    fSurf.write(np > 0);
    fFilm.write(np > 0);
    nThermoSpecies.write(np > 0);
    nFuelSpecies.write(np > 0);
    haGas.write(np > 0);
    hD.write(np > 0);
    haSurf.write(np > 0);
    haFilm.write(np > 0);
    qR.write(np > 0);    
    qD.write(np > 0);
    mFlux.write(np > 0);
    molWtSstate.write(np > 0);     
    BM.write(np > 0);

    // Write the composition fractions
    const wordList& compTypes = compModel.componentNames();

    forAll(c.linkFG(), j)
    {
        IOField<scalar> YFuel
        (
            c.newIOobject
            (
                "YFuel" + compTypes[c.linkFG()[j]],
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
            const sprayThermoParcel<ParcelType>& p = iter();             
            YFuel[i] = p.YFuel()[j];
            i++;
        }
        YFuel.write();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sprayThermoParcel<ParcelType>& p
)
{
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.TGas()
            << token::SPACE << p.TFilm()
            << token::SPACE << p.TSurf()
            << token::SPACE << p.pc()
            << token::SPACE << p.fGas()
            << token::SPACE << p.fSurf()
            << token::SPACE << p.fFilm()
            << token::SPACE << p.nThermoSpecies()
            << token::SPACE << p.nFuelSpecies()
	        << token::SPACE << p.YFuel()
            << token::SPACE << p.haGas()
            << token::SPACE << p.hD()   
            << token::SPACE << p.haSurf()
            << token::SPACE << p.haFilm()
            << token::SPACE << p.qR()
            << token::SPACE << p.qD()
            << token::SPACE << p.mFlux()
            << token::SPACE << p.molWtSstate()    
            << token::SPACE << p.BM();

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const sprayThermoParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
