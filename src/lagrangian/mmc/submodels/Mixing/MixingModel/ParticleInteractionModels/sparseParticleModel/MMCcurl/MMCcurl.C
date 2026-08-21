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

#include "MMCcurl.H"
#include "fvMesh.H"
#include "StochasticLib.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::MMCcurl<CloudType>::MMCcurl
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
Foam::MMCcurl<CloudType>::MMCcurl
(
    const MMCcurl<CloudType>& cm
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
void Foam::MMCcurl<CloudType>::mixpair
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
        scalar tauP;
        scalar tauQ;

        bool mixOn = true;

        // number of reference variables 
        const label nRefVar = pEulFields.magSqrRefVar().size();

		//Info << "nRefVar " << nRefVar << endl;

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

        if (aISO_)   //aISO mixing time scale
        {
            scalar A = 
            (
                //pEulFields.D() + pEulFields.Dt() * (dx_pq/pEulFields.DeltaE())  
				pEulFields.D() + pEulFields.Dt()
            ); 
            // denomintor for p
            scalar B = 
            (
                //qEulFields.D() + qEulFields.Dt() * (dx_pq/qEulFields.DeltaE())  
                qEulFields.D() + qEulFields.Dt()
            );

            if( A < VSMALL)
                tauP = 1e30;
            else
            {
                //tauP = CE_ * sqr(dx_pq) / (2.0 * A);
                //tauP = sqr(pEulFields.DeltaE())/(CE_*A);
		tauP = (1.0 / pEulFields.vb()) * (sqr(pEulFields.DeltaE()) / (CE_*A));
            }

            if( B < VSMALL )
                tauQ = 1e30;
            else
            {
                //tauQ = CE_ * sqr(dx_pq) / (2.0 * B );
                //tauQ = sqr(qEulFields.DeltaE()) / (CE_ * B );
		tauQ = (1.0/qEulFields.vb()) * (sqr(qEulFields.DeltaE()) / (CE_*B));
            }
        }
        else  //C&K mixing time scale
        {
			Info << "C&K called" << endl;
            const List<scalar>& magSqrfp = pEulFields.magSqrRefVar();
            const List<scalar>& magSqrfq = qEulFields.magSqrRefVar();
            forAll(magSqrfp,i)
            {
                if(magSqrfp[i]*pEulFields.DEff() < VSMALL)
                    continue;
                else
                {
                    tauPi[i] = CL_ * CE_ * beta_ *
                        sqr(p.XiR()[i] - q.XiR()[i])
                        /(magSqrfp[i]*pEulFields.DEff());

                }

                if( magSqrfq[i]*qEulFields.DEff() < VSMALL )
                    continue;
                else
                {
                    tauQi[i] = CL_ * CE_ * beta_ *
                        sqr(p.XiR()[i] - q.XiR()[i])
                        /(magSqrfq[i]*qEulFields.DEff());
                }
            }

            // mix scale from all posible mixing scales
            tauP = *(std::min_element(tauPi.begin(), tauPi.end()));

            // mix scale from all posible mixing scales
            tauQ = *(std::min_element(tauQi.begin(), tauQi.end()));
        }

        if (tauP >= 1e30 || tauQ >= 1e30)
            mixOn = false;

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
                //tauMix = max(tauP,tauQ);

            scalar mixExtent = 1.0 - exp(-deltaT / (tauMix + VSMALL));

            
//            if (!this->owner().sootingFlame())
	    if (this->owner().secondCondMixingEnabled())
            {
		// Second conditioning active: first conditioning mixes ONLY the
                // reaction progress variable phi, for every pair (no flagged/
                // unflagged distinction). Composition (Y, hA, XiC) is NOT mixed
                // here; the flagged subset mixes Y/T/hA later in second
                // conditioning. phi is reacted afterwards by W(phi) and is mixed
                // nowhere else.
                const scalar wtSum = p.wt() + q.wt();
                if (wtSum > VSMALL)
                {
                    const scalar phiAv =
                        (p.wt()*p.phi() + q.wt()*q.phi())/wtSum;
                    p.phi() += mixExtent*(phiAv - p.phi());
                    q.phi() += mixExtent*(phiAv - q.phi());
                }
	    }
	    else if (!this->owner().sootingFlame())
	    {
                particleType::mixProperties(p,q,mixExtent);
            }
            else // if it is a sooting flame
            {   
                // mix soot species at different rate
                
                // calculate mean gas diffusivity
                scalar Dgas = 0.5 * (pEulFields.D() + qEulFields.D());

                scalar Tsoot = 0.5 * (p.T() + q.T());

                // Mean viscosity
                scalar mugas = 0.5 * (pEulFields.mu() + qEulFields.mu());

                //calc Kn number Kn = rl/R_soot for all soot species >= MW200
                //rl mean free path = k*T/[sqrt(2)*pi*d^2*p]
                //d takes kinetic diameter of N2, p = 101325 [Pa]
                scalar rl = 2.3147e-10 * Tsoot;

                const wordList SpeciesNames = this->owner().composition().componentNames();

                scalarList ScaledExtent(SpeciesNames.size());

                //calc mixExtent for species
                forAll(SpeciesNames,Spi)
                {
                    if (this->owner().composition().molWt(Spi) < 200)
                        ScaledExtent[Spi] = mixExtent;
                    else //calculate Kn
                    {
                        
                        scalar sootMW = this->owner().composition().molWt(Spi);

                        scalar sootParticleMass = sootMW/1000/6.023e23; //[kg/per molecule]

                        scalar sootrho(this->owner().sootSolidDensity());

                        scalar sootV = (3.0*sootParticleMass)/(4.0*M_PI*sootrho);

                        scalar R_soot = pow(sootV,1.0/3.0);

                        scalar Kn = rl/R_soot;


                        //if Kn > a threshold, mixed as gas species (equal diffusivity)
                        if (Kn > this->owner().sootMinKn())
                        {
                            ScaledExtent[Spi] = mixExtent;
                        }
                        else
                        {   
                            // calc diffusivity of soot species
                            // D_soot = (C*K*T)/(3*pi*mu*d)
                            // C= 1 + Kn(1.257 + 0.4 exp( -0.55/(Kn/2)))
                            scalar C = 1.0 + Kn*(1.257 + 0.4*exp(-0.55/(Kn/2.0)));
                            scalar k = 1.38065e-23;

                            scalar D_soot = ((C*k*Tsoot)/(6.0*M_PI*R_soot * mugas));

                            scalar Scaled_tauMix = Dgas/D_soot * tauMix;

                            ScaledExtent[Spi] = 1.0 - exp(-deltaT / (Scaled_tauMix + VSMALL));
                        }
                    }
                }

                // ScaledExtent is a list. 
                particleType::mixProperties(p,q,mixExtent,ScaledExtent);
            } // end of if sooting flame
        } // end of if mixing
    }
}


template<class CloudType>
const Foam::scalarField Foam::MMCcurl<CloudType>::XiR0
(
    label patch,
    label patchFace,
    particle& p
)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi(); //indexes in the full set
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR(); //indexes in the reference set

    scalarField XXi(this->numXiR(),0.0);

    for (const word& varName : this->XiRNames_)
    //{
    //    if (setOfXi.Vars(XiIndexes[varName]).type()=="evolved")
    //    {
    //        XXi[XiRIndexes[varName]] = 0.;//particle.position();Need to do something when reference is shadow position
    //    }
    //    else
    //    {
    //        XXi[XiRIndexes[varName]] = setOfXi.Vars(XiIndexes[varName]).field().boundaryField()[patch][patchFace]; //setting the BCs for the reference variables on the boundary
    //    }
    //}

    {
	XXi[XiRIndexes[varName]] = setOfXi.Vars(XiIndexes[varName]).evolveMethod().initialise(patch,patchFace,p);
    }
    return XXi;

}

template<class CloudType>
const Foam::scalarField Foam::MMCcurl<CloudType>::XiR0(label celli, particle& p)
{
    const mmcVarSet& setOfXi = this->XiR();

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(this->numXiR(),0.0);

    for (const word& varName : this->XiRNames_)
    //{
    //    if (setOfXi.Vars(XiIndexes[varName]).type()=="evolved")
    //    {
    //        XXi[XiRIndexes[varName]] = 0.;//particle.position();Need to do something when reference is shadow position
    //    }
    //    else
    //    {
    //        XXi[XiRIndexes[varName]]=setOfXi.Vars(XiIndexes[varName]).field()[celli]; //assigning the cell value to the particle 
    //    }
    //}

    {
	XXi[XiRIndexes[varName]] = setOfXi.Vars(XiIndexes[varName]).evolveMethod().initialise(celli, p);
    }
    return XXi;
}


template<class CloudType>
void Foam::MMCcurl<CloudType>::printInfo()
{
    // Print info statement
    Info << "Mixing Model: " << this->modelType() << nl
         << token::TAB << "ri:      " << this->ri_<< nl
         << token::TAB << "Xi:      " << this->Xii_<< nl
         << token::TAB << "CL:      " << this->CL_<< nl
         << token::TAB << "CE:      " << this->CE_<< nl
         << token::TAB << "beta:    " << this->beta_<< nl
         << token::TAB << "aISO:    " << this->aISO_<< nl
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
}
// ************************************************************************* //
