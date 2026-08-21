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

#include "sootMMCCurl.H"
#include "fvMesh.H"
#include "StochasticLib.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::sootMMCCurl<CloudType>::sootMMCCurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    mixParticleModel<CloudType>(dict,owner, typeName, Xi),

    CL_(this->coeffDict().lookupOrDefault("CL", 0.5)),

    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),

    beta_(this->coeffDict().lookupOrDefault("beta", 3)),

    aISO_(this->coeffDict().lookupOrDefault("aISO",true)),

    meanTimeScale_(this->coeffDict().lookup("meanTimeScale"))
{
    printInfo();
}


template <class CloudType>
Foam::sootMMCCurl<CloudType>::sootMMCCurl
(
    const sootMMCCurl<CloudType>& cm
)
:
    mixParticleModel<CloudType>(cm),

    CL_(this->coeffDict().lookupOrDefault("CL", 0.5)),

    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),

    beta_(this->coeffDict().lookupOrDefault("beta", 3)),

    aISO_(this->coeffDict().lookupOrDefault("aISO",true))
{
    printInfo();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::sootMMCCurl<CloudType>::mixpair
(
    particleType& p,
    const eulerianFieldData& pEulFields,
    particleType& q,
    const eulerianFieldData& qEulFields,
    scalar& deltaT
)
{

    //- If combined weights of p and q > 0
    if(p.wt() + q.wt() > 0)
    {
        //- Mixing time scale
        scalar tauP;
        scalar tauQ;
    
        scalar tauPsoot;
        scalar tauQsoot;

        bool mixOn = true;

        // number of reference variables 
        const label nRefVar = pEulFields.magSqrRefVar().size();

        List<scalar> tauPi(nRefVar,1e30);

        List<scalar> tauQi(nRefVar,1e30);

        // calculate the distance in physical space
        scalar dx_pq = 
            sqrt
            (
                sqr(pEulFields.position().x()-qEulFields.position().x())
              + sqr(pEulFields.position().y()-qEulFields.position().y())
              + sqr(pEulFields.position().z()-qEulFields.position().z())
            );
            
        p.dx() = dx_pq;
        q.dx() = dx_pq;

        forAll(p.dXiR(),i)
        {
            // Set distance in each reference space
            p.dXiR()[i] = mag(p.XiR()[i] - q.XiR()[i]);
            q.dXiR()[i] = mag(p.XiR()[i] - q.XiR()[i]);
        }

        //- Compute mixing time scales for particle p
        scalar A = 
        (
            pEulFields.D() + pEulFields.Dt() * (dx_pq/pEulFields.DeltaE())  
        ); 
        
        scalar B = 
        (
            qEulFields.D() + qEulFields.Dt() * (dx_pq/qEulFields.DeltaE())  
        );

        scalar A4soot = 
        (
            pEulFields.Dt() * (dx_pq/pEulFields.DeltaE())  
        ); 
        
        scalar B4soot = 
        (
            qEulFields.Dt() * (dx_pq/qEulFields.DeltaE())  
        );

        if( A < VSMALL)
            tauP = 1e30;
        else
        {
            tauP = CE_ * sqr(dx_pq) / (2.0 * A);
        }

        if( B < VSMALL )
            tauQ = 1e30;
        else
        {
            tauQ = CE_ * sqr(dx_pq) / (2.0 * B );
        }

        if( A4soot < VSMALL)
            tauPsoot = 1e30;
        else
        {
            tauPsoot = CE_ * sqr(dx_pq) / (2.0 * A4soot);
        }

        if( B4soot < VSMALL )
            tauQsoot = 1e30;
        else
        {
            tauQsoot = CE_ * sqr(dx_pq) / (2.0 * B4soot );
        }


        if(tauP >= 1e30 || tauQ >= 1e30)
            mixOn = false;

        if (mixOn && p.wt()+q.wt() > 1e-200)
        {
            scalar tauMix = 0.0;
            scalar tauMixSoot = 0.0;

            if (meanTimeScale_)
            {
                tauMix = 2.0
                   /(
                        1.0/(tauP + VSMALL) 
                      + 1.0/(tauQ + VSMALL)
                    );
                tauMixSoot = 2.0
                   /(
                        1.0/(tauPsoot + VSMALL) 
                      + 1.0/(tauQsoot + VSMALL)
                    );
            }
            else
            {
                tauMix = min(tauP,tauQ);
                tauMixSoot = min(tauPsoot, tauQsoot);
            }

            scalar mixExtent = 1.0 - exp(-deltaT / (tauMix + VSMALL));
            scalar mixExtentSoot = 1.0 - exp(-deltaT / (tauMixSoot + VSMALL));

            
            if(this->owner().soot().turbMixing())
                particleType::mixProperties(p,q,mixExtent,mixExtentSoot);
            else 
                particleType::mixProperties(p,q,mixExtent);
        }
    }
}


template<class CloudType>
const Foam::scalarField Foam::sootMMCCurl<CloudType>::XiR0
(
    label patch,
    label patchFace
)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(),0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIndexes[varName]).type()=="evolved")
        {
            XXi[XiRIndexes[varName]] = 0.;//particle.position();Need to do something when reference is shadow position
        }
        else
        {
            XXi[XiRIndexes[varName]] = setOfXi.Vars(XiIndexes[varName]).field().boundaryField()[patch][patchFace];
        }
    }
    return XXi;

}

template<class CloudType>
const Foam::scalarField Foam::sootMMCCurl<CloudType>::XiR0(label celli)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(),0.0);

    for (const word& varName : this->XiRNames_)
    {
        if (setOfXi.Vars(XiIndexes[varName]).type()=="evolved")
        {
            XXi[XiRIndexes[varName]] = 0.;//particle.position();Need to do something when reference is shadow position
        }
        else
        {
            XXi[XiRIndexes[varName]]=setOfXi.Vars(XiIndexes[varName]).field()[celli];
        }
    }
    return XXi;
}


template<class CloudType>
void Foam::sootMMCCurl<CloudType>::printInfo()
{
    // Print info statement
    Info << "Mixing Model: " << this->modelType() << nl
         << token::TAB << "ri:      " << this->ri_<< nl
         << token::TAB << "Xi:      " << this->Xii_<< nl
         << token::TAB << "CL:      " << this->CL_<< nl
         << token::TAB << "CE:      " << this->CE_<< nl
         << token::TAB << "beta:    " << this->beta_<< nl
         << token::TAB << "aISO:    " << Switch(true)<< nl
         << token::TAB << "---------------------------"<< endl;
    Info << token::TAB << "General Mixing Rules:    " << nl
         << token::TAB << "Particle pairing method: "  
         << this->pairingMethod_ << endl;
    if(meanTimeScale_)
        Info << token::TAB 
             << "Use harmonic mean of two particles' mixing time scales" 
             << endl;
    else
        Info << token::TAB
             << "Use the min of two particles' mixing time scales" 
             << endl;
    
    if (!aISO_)
        Info << token::TAB
             << "WARNING: Currently only aISO model is supported. Will use "
             << "aISO model!"
             << endl;
}
// ************************************************************************* //
