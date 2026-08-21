/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "Curl.H"
#include "fvMesh.H"
#include "StochasticLib.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::Curl<CloudType>::Curl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    denseMixParticleModel<CloudType>(dict,owner, typeName),
    
    CL_(this->coeffDict().template lookupOrDefault<scalar>("CL", 0.5)),
    
    CE_(this->coeffDict().template lookupOrDefault<scalar>("CE", 0.1)),

    XiR_(Xi)
    
{
    // Print the dictionary
    Info << "Mixing dictionary "<<this->coeffDict()<<endl;  
}


template <class CloudType>
Foam::Curl<CloudType>::Curl
(
    const Curl<CloudType>& cm
)
:
    denseMixParticleModel<CloudType>(cm),
    
    CL_(this->coeffDict().template lookupOrDefault<scalar>("CL", 0.5)),
    
    CE_(this->coeffDict().template lookupOrDefault<scalar>("CE", 0.1)),

    XiR_(cm.XiR_)
    
{
    // Print the dictionary
    Info << "Mixing dictionary "<<this->coeffDict()<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::Curl<CloudType>::~Curl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::Curl<CloudType>::mixpair
(
    particleType& p, 
    const eulerianFieldParticleValues& pEulField,
    particleType& q,         
    const eulerianFieldParticleValues& qEulField,
    scalar& deltaT
)
{
    //- If combined weights of p and q > 0
    if(p.wt() + q.wt() > 0)
    {
        //- Mixing time scale
        scalar tauP;
      
        scalar tauQ;
        
        bool mixOn = true;

        if(pEulField.DEff() < VSMALL)
        {
            tauP = 1e30;
            
            mixOn = false;
        }
        else
        {
            tauP = CL_ * CE_ * pEulField.sqrFilter() / pEulField.DEff();
        }

        if(qEulField.DEff() < VSMALL)
        {   
            tauQ = 1e30;
            
            mixOn = false;
        }
        else
        {
            tauQ = CL_ * CE_ * qEulField.sqrFilter() / qEulField.DEff();
        }

        if (mixOn)
        {
            scalar tauMix = max(tauP,tauQ);
      
            scalar mixExtent = 1 - exp(-deltaT / (tauMix + VSMALL));         
            
            // Function cascade to mix scalars
            mixScalarProperties(p,q,mixExtent);
        }
    }
}


template<class CloudType>
void Foam::Curl<CloudType>::SmixList
(
    DynamicList<particleType*>& particleList,
    DynamicList<eulerianFieldParticleValues>& eulFields
)
{
    scalar deltaT = this->owner().mesh().time().deltaT().value();
      
    if(particleList.size() > 3) 
    {
        // Shuffle list to enable random mixing of particles 
        // As both lists, eulerianFieldsList_ and particlePtrList have to 
        // shuffeled in the same way an index list is used
        List<label> indices(particleList.size());
        std::iota(indices.begin(),indices.end(),0);
        std::random_shuffle(indices.begin(),indices.end());
        
        for(int i=0; i<std::floor(0.5*particleList.size()); i++)
        {
            label p = indices[i*2];
                
            label q = indices[i*2+1];
                
            mixpair(*particleList[p],eulFields[p],*particleList[q],eulFields[q],deltaT);
        }
    }
}

template<class CloudType>
inline const Foam::scalarField Foam::Curl<CloudType>::XiR0(label patch, label patchFace) 
{  
    return scalarField();
}

template<class CloudType>
inline const Foam::scalarField Foam::Curl<CloudType>::XiR0(label celli) 
{
    return scalarField();
}


template<class CloudType>
void Foam::Curl<CloudType>::mixScalarProperties
(
    particleType& p,
    particleType& q,
    scalar mixExtent
)
{
    typedef typename CloudType::particleType particleType;
    
    particleType::mixProperties(p,q,mixExtent);
}


//- return mmcVarSet
template<class CloudType>
const Foam::mmcVarSet& Foam::Curl<CloudType>::XiR() const
{
    return XiR_;
}

// ************************************************************************* //
