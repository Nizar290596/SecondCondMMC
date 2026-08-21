/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DropletSprayThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::DropletSprayThermoParcel<ParcelType>::propertyList_ =
    Foam::DropletSprayThermoParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::DropletSprayThermoParcel<ParcelType>::sizeofFields
(
    sizeof(DropletSprayThermoParcel<ParcelType>)
  - offsetof(DropletSprayThermoParcel<ParcelType>, T_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DropletSprayThermoParcel<ParcelType>::DropletSprayThermoParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    T_(0.0),
    Cp_(0.0),
    jFuel_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> T_ >> Cp_ >> mass0_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &T_);
            readRawScalar(is, &Cp_);
            readRawScalar(is, &mass0_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&T_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);

    initStatisticalSampling();
}


template<class ParcelType>
template<class CloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::readFields(CloudType& c)
{
    const bool readOnProc = c.size();

    ParcelType::readFields(c);

    IOField<scalar> T(c.newIOobject("T", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, T);

    IOField<scalar> Cp(c.newIOobject("Cp", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, Cp);

    IOField<scalar> mass0(c.newIOobject("mass0", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, mass0);


    label i = 0;
    for (DropletSprayThermoParcel<ParcelType>& p : c)
    {
        p.T_ = T[i];
        p.Cp_ = Cp[i];
        p.mass0_ = mass0[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<scalar> T(c.newIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.newIOobject("Cp", IOobject::NO_READ), np);
    IOField<scalar> mass0(c.newIOobject("mass0", IOobject::NO_READ), np);

    label i = 0;
    for (const DropletSprayThermoParcel<ParcelType>& p : c)
    {
        T[i] = p.T_;
        Cp[i] = p.Cp_;
        mass0[i] = p.mass0_;

        ++i;
    }

    T.write(writeOnProc);
    Cp.write(writeOnProc);
    mass0.write(writeOnProc);
}


template<class ParcelType>
void Foam::DropletSprayThermoParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("T", T_);
    writeProp("Cp", Cp_);
    writeProp("mass", mass0_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readFields(c);

    if (!c.size()) return;

    auto& T = cloud::lookupIOField<scalar>("T", obr);
    auto& Cp = cloud::lookupIOField<scalar>("Cp", obr);
    auto& mass0 = cloud::lookupIOField<scalar>("mass0", obr);

    label i = 0;
    for (DropletSprayThermoParcel<ParcelType>& p : c)
    {
        p.T_ = T[i];
        p.Cp_ = Cp[i];
        p.mass0_ = mass0[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& T = cloud::createIOField<scalar>("T", np, obr);
    auto& Cp = cloud::createIOField<scalar>("Cp", np, obr);
    auto& mass0 = cloud::createIOField<scalar>("mass", np, obr);

    label i = 0;
    for (const DropletSprayThermoParcel<ParcelType>& p : c)
    {
        T[i] = p.T_;
        Cp[i] = p.Cp_;
        mass0[i] = p.mass0_;

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DropletSprayThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.Cp()
            << token::SPACE << p.mass0();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.T_),
            DropletSprayThermoParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
