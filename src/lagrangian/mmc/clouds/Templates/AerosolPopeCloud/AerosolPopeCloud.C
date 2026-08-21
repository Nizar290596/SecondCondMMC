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

#include "AerosolPopeCloud.H"
#include "SynthesisModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::setModels()
{
    synthesisModel_.reset
    (
        SynthesisModel<AerosolPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::cloudReset(AerosolPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    synthesisModel_.reset(c.synthesisModel_ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AerosolPopeCloud<CloudType>::AerosolPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        U,
        DEff,
        rho,
        gradRho,
        Xi,
        false,
        true
    ),

    aerosolPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    synthesisModel_(nullptr)
{
    Info << nl << "Creating aerosol Pope Particle Cloud." << nl << endl;
    
    setModels();
    
    Info << nl << "Aerosol model constructed." << endl;
    
    setEulerianStatistics();    
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Aerosol Pope particle cloud data from file." << endl;
        
            particleType::aerosolParticleIOType::readFields(*this,this->synthesis());
        }
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of aerosol Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
            }
        }
    }
}



template<class CloudType>
Foam::AerosolPopeCloud<CloudType>::AerosolPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        mesh,
        Xi,
        false,          // Needs to be false as top level calls the initAtCnstr
        readFields
    ),

    aerosolPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    synthesisModel_(nullptr)
{
    Info << nl << "Creating aerosol Pope Particle Cloud." << nl << endl;
    
    setModels();
    
    Info << nl << "Aerosol model constructed." << endl;
    
    setEulerianStatistics();    
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Aerosol Pope particle cloud data from file." << endl;
        
            particleType::aerosolParticleIOType::readFields(*this,this->synthesis());
        }
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial release of aerosol Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AerosolPopeCloud<CloudType>::~AerosolPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle, mass, wt, patchI, patchFace,iniRls);
    

    // Set static properties if they are uninitialized
    if 
    (
        particleType::physicalAerosolPropertyNames_.empty()
     && particleType::psdPropertyNames_.empty()
    )
    {
        particleType::physicalAerosolPropertyNames_ = 
            synthesis().physicalAerosolPropertyNames();
        particleType::psdPropertyNames_ = 
            synthesis().psdPropertyNames();
    }
    
    particle.physicalAerosolProperties().setSize(synthesis().physicalAerosolPropertyNames().size(), 0.0);
    
    particle.psdProperties().setSize(synthesis().psdPropertyNames().size(), 0.0);

    particle.initStatisticalSampling();
}


template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::setEulerianStatistics()
{
    const wordList& physicalAerosolPropertyNames = synthesis().physicalAerosolPropertyNames();
    
    forAll(physicalAerosolPropertyNames, physicalAerosolProperty)
    {
        if (this->eulerianStatsDict().found(physicalAerosolPropertyNames[physicalAerosolProperty]))
        {
            const dimensionSet dim = dimless;
        
            this->eulerianStats().newProperty(physicalAerosolPropertyNames[physicalAerosolProperty],dim);
        } 
    }
    
    const wordList& psdPropertyNames = synthesis().psdPropertyNames();
    
    forAll(psdPropertyNames, psdProperty)
    {
        if (this->eulerianStatsDict().found(psdPropertyNames[psdProperty]))
        {
            const dimensionSet dim = dimless;
        
            this->eulerianStats().newProperty(psdPropertyNames[psdProperty],dim);
        } 
    }
}


template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();
    
    const wordList& physicalAerosolPropertyNames = synthesis().physicalAerosolPropertyNames();

    const wordList& psdPropertyNames = synthesis().psdPropertyNames();

    forAllIters(*this, iter)
    {
        this->eulerianStats().findCell(iter().position());

        forAll(physicalAerosolPropertyNames, physicalAerosolProperty)
        {
            if (this->eulerianStatsDict().found(physicalAerosolPropertyNames[physicalAerosolProperty]))
                this->eulerianStats().calculate(physicalAerosolPropertyNames[physicalAerosolProperty],iter().wt(),iter().physicalAerosolProperties()[physicalAerosolProperty]);
        }

        forAll(psdPropertyNames, psdProperty)
        {
            if (this->eulerianStatsDict().found(psdPropertyNames[psdProperty]))
                this->eulerianStats().calculate(psdPropertyNames[psdProperty],iter().wt(),iter().psdProperties()[psdProperty]);
        }
    }
}


template<class CloudType>
void Foam::AerosolPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::aerosolParticleIOType::writeFields(*this, this->synthesis());
    }
}


// ************************************************************************* //
