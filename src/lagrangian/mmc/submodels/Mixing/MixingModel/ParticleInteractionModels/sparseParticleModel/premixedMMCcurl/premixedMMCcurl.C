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

#include "premixedMMCcurl.H"
#include "fvMesh.H"
#include "StochasticLib.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::premixedMMCcurl<CloudType>::premixedMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    premixedMixParticleModel<CloudType>(dict,owner, typeName, Xi),
    
    mixTimeScale_(this->coeffDict().lookup("mixTimeScale")),

    CL_(this->coeffDict().lookupOrDefault("CL", 0.1)),

    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),
     
    Ka_(readScalar(this->coeffDict().lookup("KarlovitzNum"))), 
    
    alpha_(this->coeffDict().lookupOrDefault("alpha", 1.0)), 

    limitkdTree_(this->coeffDict().lookupOrDefault("limitkdTree", true)),
    
    meanTimeScale_(this->coeffDict().lookup("meanTimeScale"))
{
    printInfo();
}


template <class CloudType>
Foam::premixedMMCcurl<CloudType>::premixedMMCcurl
(
    const premixedMMCcurl<CloudType>& cm
)
:
    premixedMixParticleModel<CloudType>(cm),
    
    mixTimeScale_(this->coeffDict().lookup("mixTimeScale")),

    CL_(this->coeffDict().lookupOrDefault("CL", 0.1)),

    CE_(this->coeffDict().lookupOrDefault("CE", 0.1)),
      
    Ka_(readScalar(this->coeffDict().lookup("KarlovitzNum"))),
    
    alpha_(this->coeffDict().lookupOrDefault("alpha", 1.0)), 

    limitkdTree_(this->coeffDict().lookupOrDefault("limitkdTree", true))   
{
    printInfo();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::premixedMMCcurl<CloudType>::mixpair
(
    particleType& p,
    const eulerianFieldData& pEulFields,
    particleType& q,
    const eulerianFieldData& qEulFields,
    scalar& deltaT
)
{

    //- If combined weights of p and q > 0
    if (p.wt() + q.wt() > 0)
    {
        //- Mixing time scale
        scalar tauP(1e30);
        scalar tauQ(1e30);

        bool mixOn = true;

        // number of reference variables 
        const label nRefVar = pEulFields.magSqrRefVar().size();

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

        //- Compute mixing time scales for particles p and q
 
        if (mixTimeScale_ == "LPF_aISO")   // LPF_aISO blended mixing time scale
        {
        
            scalar tauPlam(1e30);
            scalar tauQlam(1e30);
            scalar tauPturb(1e30);
            scalar tauQturb(1e30);

            // - LPF mixing time scale (mixTime values come from solver)
            tauPlam = CL_ * p.mixTime();
            tauQlam = CL_ * q.mixTime();            
            
            // - Compute aISO mixing time scale 
            scalar A = 
            (
                pEulFields.D() + pEulFields.Dt() * (dx_pq/pEulFields.DeltaE())  
            ); 
            
            scalar B = 
            (
                qEulFields.D() + qEulFields.Dt() * (dx_pq/qEulFields.DeltaE())  
            );

            if ( A < VSMALL )
                tauPturb = 1e30;
            else
            {
                tauPturb = CE_ * sqr(dx_pq) / (2.0 * A);
            }

            if ( B < VSMALL )
                tauQturb = 1e30;
            else
            {
                tauQturb = CE_ * sqr(dx_pq) / (2.0 * B);
            }
            
            // - Blending function for the laminar and turbulent mixing time scales
            
            scalar blendFuncP(0.0);
            scalar blendFuncQ(0.0);
            scalar cRLXp(0.0);
            scalar cRLXq(0.0);
            
            if (nRefVar > 1) 
            {
                FatalError << "Mixing time scale LPF_aISO is implemented only for a single reference variable." 
                           << exit(FatalError);
            }
            else 
            {
                cRLXp = p.XiR()[0];
                cRLXq = q.XiR()[0];
            }
            
            if (Ka_ > 1.0)
            {
                blendFuncP = 2.0 * sqrt(cRLXp * (1.0 - cRLXp)) / (pow(Ka_, alpha_));
                blendFuncQ = 2.0 * sqrt(cRLXq * (1.0 - cRLXq)) / (pow(Ka_, alpha_));
            }
            else 
            {
                blendFuncP = 2.0 * sqrt(cRLXp * (1.0 - cRLXp));
                blendFuncQ = 2.0 * sqrt(cRLXq * (1.0 - cRLXq));
            }
            
            if (tauPturb >= 1e30 || tauPlam >= 1e30)
                tauP = 1e30;
            else if ((tauPlam <= tauPturb) && (blendFuncP > 0.05))
                tauP = tauPlam;
            else
            {
                tauP = 1.0 / ( blendFuncP / (tauPlam + VSMALL) + (1.0 - blendFuncP) / (tauPturb + VSMALL));
            }

            if (tauQturb >= 1e30 || tauQlam >= 1e30)
                tauQ = 1e30;
            else if ((tauQlam <= tauQturb) && (blendFuncQ > 0.05))
                tauQ = tauQlam;
            else
            {
                tauQ = 1.0 / ( blendFuncQ / (tauQlam + VSMALL) + (1.0 - blendFuncQ) / (tauQturb + VSMALL));
            }
            
        }
        else if (mixTimeScale_ == "LPF")   // mixing time scale for Laminar Premixed Flames (mixTime values come from solver)
        {
            tauP = CL_ * p.mixTime();
            tauQ = CL_ * q.mixTime();
        }
        else 
        {
            FatalError << "Selected mixing time scale model does not exist."<<nl
                       << "Valid models are: LPF_aISO, LPF"
                       << exit(FatalError);
        } 

        if (tauP >= 1e30 || tauQ >= 1e30)
            mixOn = false;
            
        // - Do not mix particles that are far in composition space (dXiR > 2*Ximi)
        if (limitkdTree_) 
        {
            forAll(p.dXiR(),i)
            {
                if (p.dXiR()[i] > 2.0 * this->Xii_[i])
                {
                    mixOn = false;
                    p.mixTime() = -1.0; // - Here values are set to -1 in order to count percentage of not mixing particles in postprocessing 
                    q.mixTime() = -1.0;
                    p.mixExt() = -1.0;  
                    q.mixExt() = -1.0;  
                }
            }
        }

        if (mixOn && p.wt()+q.wt() > 1e-200)
        {
            scalar tauMix = 0.0;

            if (meanTimeScale_)
                tauMix = 2.0
                   /(
                        1.0/(tauP + VSMALL) 
                      + 1.0/(tauQ + VSMALL)
                    );
            else
                tauMix = min(tauP,tauQ);

            scalar mixExtent = 1.0 - exp(-deltaT / (tauMix + VSMALL));
            
            particleType::mixProperties(p,q,mixExtent);
            
            p.mixTime() = tauMix; 
            q.mixTime() = tauMix;  
            p.mixExt() = mixExtent; 
            q.mixExt() = mixExtent;
        } // end of if mixing
    }
}


template<class CloudType>
const Foam::scalarField Foam::premixedMMCcurl<CloudType>::XiR0
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
const Foam::scalarField Foam::premixedMMCcurl<CloudType>::XiR0(label celli)
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
void Foam::premixedMMCcurl<CloudType>::printInfo()
{
    // Print info statement
    
    if (mixTimeScale_ == "LPF_aISO") 
    {
        Info << "Mixing Time Scale: LPF_aISO (turbulent premixed flame)" << nl
             << token::TAB << "CL:                " << this->CL_<< nl
             << token::TAB << "CE:                " << this->CE_<< nl
             << token::TAB << "Karlovitz number:  " << this->Ka_<< nl             
             << token::TAB << "alpha:             " << this->alpha_<< endl;
    }
    else if (mixTimeScale_ == "LPF") 
    {
        Info << "Mixing Time Scale: LPF (laminar premixed flame)" << nl
             << token::TAB << "CL:                " << this->CL_<< endl;
    }
    else
    {
        FatalError << "Selected mixing time scale model does not exist."<< nl
                   << "Valid models are: LPF_aISO, LPF"
                   << exit(FatalError);
    }
    Info << token::TAB << "ri:                " << this->ri_<< nl
         << token::TAB << "Xi:                " << this->Xii_<< nl
         << token::TAB << "---------------------------"<< endl;
    Info << token::TAB << "General Mixing Rules:    " << nl
         << token::TAB << "Particle pairing method: "  
         << this->pairingMethod_ << endl;
    if (meanTimeScale_)
        Info << token::TAB 
             << "Use harmonic mean of two particles' mixing time scales" 
             << endl;
    else
        Info << token::TAB
             << "Use the min of two particles' mixing time scales" 
             << endl;
    if (this->ATFkd_)
        Info << token::TAB 
             << "ATF is accounted for in kdTree pairing" 
             << endl;
    else
        Info << token::TAB
             << "ATF is not accounted for in kdTree pairing"
             << endl;
             
    const Switch ATFmixTimeOn_(this->coeffDict().lookup("ATFmixTime"));
    if (ATFmixTimeOn_)
        Info << token::TAB 
             << "ATF is accounted for in LPF mixing time scale distribution" 
             << endl;
    else
        Info << token::TAB
             << "ATF is not accounted for in LPF mixing time scale distribution"
             << endl;
    
    if (limitkdTree_)
        Info << token::TAB 
             << "Particles with dXiR > 2*Ximi are NOT mixed" 
             << endl;
    else
        Info << token::TAB
             << "No limiting of mixing based on dXiR" 
             << endl;
}
// ************************************************************************* //
