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

#include "AerosolPopeParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::AerosolPopeParticle<ParticleType>::propertyList_ =
    Foam::AerosolPopeParticle<ParticleType>::propertyList();

template<class ParticleType>
Foam::List<Foam::word> 
Foam::AerosolPopeParticle<ParticleType>::physicalAerosolPropertyNames_= 
    Foam::List<Foam::word>(0);
    
template<class ParticleType>
Foam::List<Foam::word> 
Foam::AerosolPopeParticle<ParticleType>::psdPropertyNames_= 
    Foam::List<Foam::word>(0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::AerosolPopeParticle<ParticleType>::AerosolPopeParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParticleType(mesh, is, readFields, newFormat),
    physicalAerosolProperties_(0),
    psdProperties_(0)
{
    if (readFields)
    {
        DynamicList<scalar> physicalAerosolProperties;
        
        is  >> physicalAerosolProperties;
        
        physicalAerosolProperties_.transfer(physicalAerosolProperties);

        DynamicList<scalar> psdProperty;
        
        is  >> psdProperty;
        
        psdProperties_.transfer(psdProperty);
    }

    initStatisticalSampling();

    // Check state of Istream
    is.check
    (
        "AerosolPopeParticle::AerosolPopeParticle(const polyMesh&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * * Read / Write * * * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::AerosolPopeParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);
}


template<class ParticleType>
template<class CloudType, class AerosolType>
void Foam::AerosolPopeParticle<ParticleType>::readFields
(
    CloudType& c,
    const AerosolType& aerosolModel
)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);

    // Get names and sizes for each physical aerosol and PSD property
    const wordList& physicalAerosolPropertyNames = aerosolModel.physicalAerosolPropertyNames();
    const wordList& propertyNames = aerosolModel.psdPropertyNames();    
    
    const label nPhysicalAerosolProperties = physicalAerosolPropertyNames.size();
    const label nProperties = propertyNames.size();
    
    
    // Set storage for each PSD properties for each particle
    forAllIters(c, iter)
    {
        AerosolPopeParticle<ParticleType>& p = iter();
        p.physicalAerosolProperties_.setSize(nPhysicalAerosolProperties, 0.0);
        p.psdProperties_.setSize(nProperties, 0.0);
    }
    
    // Populate physical aerosol properties for each particle
    forAll(physicalAerosolPropertyNames, j)
    {
        IOField<scalar> physicalAerosolProperties
        (
            c.newIOobject
            (
                physicalAerosolPropertyNames[j],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            AerosolPopeParticle<ParticleType>& p = iter();
            p.physicalAerosolProperties_[j] = physicalAerosolProperties[i];
            
            i++;
        }
    }
    
    // Populate PSD properties for each particle
    forAll(propertyNames, j)
    {
        IOField<scalar> psdProperties
        (
            c.newIOobject
            (
                propertyNames[j],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIters(c, iter)
        {
            AerosolPopeParticle<ParticleType>& p = iter();
            p.psdProperties_[j] = psdProperties[i];
            
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
void Foam::AerosolPopeParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);
}


template<class ParticleType>
template<class CloudType, class AerosolType>
void Foam::AerosolPopeParticle<ParticleType>::writeFields
(
    const CloudType& c,
    const AerosolType& aerosolModel
)
{
    ParticleType::writeFields(c);

    label np = c.size();
    if (np > 0)
    {
        // Write the physical aerosol and PSD property
        const wordList& physicalAerosolPropertyNames = aerosolModel.physicalAerosolPropertyNames();
        
        forAll(physicalAerosolPropertyNames, j)
        {
            IOField<scalar> physicalAerosolProperties
            (
                c.newIOobject
                (
                    physicalAerosolPropertyNames[j],
                    IOobject::NO_READ
                ),
                np
            );
      
            label i = 0;
            forAllConstIter
            (
                typename Cloud<AerosolPopeParticle<ParticleType> >,
                c,
                iter
            )
            {
                const AerosolPopeParticle<ParticleType>& p = iter();
                
                physicalAerosolProperties[i] = p.physicalAerosolProperties()[j];

                i++;
            }
            
            physicalAerosolProperties.write();
        }        
        
        const wordList& propertyNames = aerosolModel.psdPropertyNames();
        
        forAll(propertyNames, j)
        {
            IOField<scalar> psdProperties
            (
                c.newIOobject
                (
                    propertyNames[j],
                    IOobject::NO_READ
                ),
                np
            );
      
            label i = 0;
            forAllConstIter
            (
                typename Cloud<AerosolPopeParticle<ParticleType> >,
                c,
                iter
            )
            {
                const AerosolPopeParticle<ParticleType>& p = iter();
                
                psdProperties[i] = p.psdProperties()[j];

                i++;
            }
            
            psdProperties.write();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const AerosolPopeParticle<ParticleType>& p
)
{
    os  << static_cast<const ParticleType&>(p)
        << token::SPACE << p.physicalAerosolProperties()
        << token::SPACE << p.psdProperties();
        
    // Check state of Ostream
    os.check
    (
       "Ostream& operator<<(Ostream&, const AerosolPopeParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //
