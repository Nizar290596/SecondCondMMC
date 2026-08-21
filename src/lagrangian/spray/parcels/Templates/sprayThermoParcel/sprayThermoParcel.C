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
#include "physicoChemicalConstants.H"

#include "interpolationLookUpTable.H"
#include "basicAerosolReactingPopeCloud.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    tetIndices tetIs = this->currentTetIndices();

    pc_ = td.PInterp().interpolate(this->coordinates(), tetIs);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const label celli = this->cell();
    td.Uc() += cloud.UTrans()[celli]/this->massCell(td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::calcSurfaceProperties
(
        TrackCloudType& cloud,
        trackingData& td
) 
{
    const scalar LH = cloud.constProps().LH(); 

    //link to gaseous species ID for the liquid species
    scalar linkFGSize = cloud.linkFG().size();

    //-Surface Mixture fraction   
    molWtSstate_ = cloud.constProps().molWtSstate0();      
    label iterate = 1;  
    scalar epsilon = 0.00125;

    while(iterate == 1)
    {
        const scalar molWtSstate_old =  molWtSstate_;
        noMassTransfer_ = 0;

        scalarField YFuelSurf;
        YFuelSurf.setSize(linkFGSize,0.0);

        scalar YFuelSurfSum = 0.0;

        YFuelSurf = cloud.evaporation().YFuelVap(YFuel_,TD_, pc_, LH, molWtSstate_);      //- YFuelSurf: scalarField with YSurf for each fuel component
        forAll(YFuelSurf,nfs) YFuelSurfSum += YFuelSurf[nfs];                             //  (sum(YFuelSurf != 0), since it'S related to the gas composition but of size of the liquid composition)
        YFuelVap_ = YFuelSurf / YFuelSurfSum;                                             //- composition of the fuel vapor, this is equal to YFuel_ for single component fuels
                                                                                          //  or multi component fuels with no diff evap

        fSurf_ = cloud.envelope().fSurface(YGas_,YFuelVap_,fGas_,YFuelSurf);
        
        if(fSurf_ <= fGas_) //- if evaporation results in YSurf_Fuel < YGas_Fuel set fSurf to fGas and disable mass transfer
        {                   //  (not good: the treatment of fSurf < fGas needs to be improved)
            fSurf_ = fGas_ + 0.001;
            noMassTransfer_ = 1;
        }

        //return composition state at surface (between YFuelVap and YGas)
        YSurf_ = cloud.envelope().YSurface(YGas_,YFuelVap_,fSurf_);

        molWtSstate_ = cloud.pope().composition().mixtureMW(YSurf_) ;

        if(fabs((molWtSstate_old - molWtSstate_)/molWtSstate_old) < epsilon) iterate = 0;
    }

    this->particleH(cloud,td);

    haSurf_ = cloud.envelope().haSurface(YSurf_,pc_,TD_,haGas_,haFuelVap_,fSurf_);

    TSurf_ = cloud.envelope().TSurface(YSurf_,pc_,TD_,TGas_,fSurf_,haSurf_);
}

//- Calculate hD_
template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::particleH
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    haFuelVap_ = cloud.pope().composition().particleMixture(this->YLiqToYGas(cloud,td,YFuelVap_)).Ha(pc_,TD_); // enthalpy of the fuels vapour
    
    scalar haFuel = cloud.pope().composition().particleMixture(this->YLiqToYGas(cloud,td,YFuel_)).Ha(pc_,TD_);
    hD_ = haFuel - cloud.constProps().LH(); // enthalpy of the liquid fuel
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::calcThermophysicalValues
(
    TrackCloudType& cloud,
    trackingData& td
)  
{
    //- Filmstate mixture fraction
    fFilm_ = (2.0*fSurf_ + fGas_)/3.0;
    YFilm_ = cloud.envelope().YFilm(YGas_,YFuelVap_,fFilm_);    

    //*********************Film state properties*********************
    //- Enthalpy of film state: independent of envelope
    haFilm_ = (2.0*haSurf_ + haGas_)/3.0;   //- According to two-third rulei
    TFilm_ = cloud.pope().composition().particleMixture(YFilm_).THa(haFilm_, pc_, 2*TGas_); //2*TGas: prevents T-range errors during iteration

    //-Calculate transport properties********************************  
    //-Calculate Filmstate density and viscosity

    rhoFilm_ = cloud.pope().composition().particleMixture(YFilm_).rho(pc_, TFilm_);

    const scalar As = 1.67212E-06;
    const scalar Ts = 170.672; 
    muFilm_ = As * sqrt(TFilm_) / (1.0 + Ts / TFilm_);// Sutherland's model for calculating film state viscosity

    //Filmstate Reynolds number
    const scalar ReFilm = this->Re(rhoFilm_, this->U(), td.Uc(), this->d(), muFilm_);   

    //- Schmidt number from constantProperties
    const scalar  ScFilm = cloud.constProps().Sc();;

    //-Lewis number form constantProperties
    const scalar LeFilm = cloud.constProps().Le();

    //- Mass diffusivity (D)
    const scalar mDiffFilm = muFilm_/(rhoFilm_*ScFilm); 

    //- Thermal diffusivity (alpha)
    const scalar thDiffFilm = mDiffFilm * LeFilm;  
            
    // Prandtl number 
    const scalar PrFilm = this->PrFilm(muFilm_, rhoFilm_, thDiffFilm);

    // Nusselt number  
    //scalar NuFilm = this->NuFilm(ReFilm, PrFilm);  
    const scalar NuFilm = 2.0 + 0.6 * sqrt(ReFilm) * pow(PrFilm, (1.0 / 3.0));

    // Sherwood number
    //const scalar ShFilm = this->ShFilm(ReFilm, Sc);
    const scalar ShFilm = 2.0 + 0.6 * sqrt(ReFilm) * pow(ScFilm, (1.0 / 3.0));
  
    //- Constant Z = Le^-1 * Sh / Nu.   NOTE: Z = 1 and BM=BH for unity Lewis No.
    ZFilm_ = ShFilm / (NuFilm * LeFilm);    
        
    // Spalding transfer number for Mass
    BM_ = this->BM(fGas_, fSurf_); 

    if (BM() == -1)
    {
        mFlux_ = 0.0; 
        qD_    = 0.0 ;
        qR_    = 0.0;
    }
    else 
    {  
        //- Sherwood number correction factor
        const scalar FM = pow((1 + BM()), 0.7) * log(1.0 + BM()) / BM();

        //- Corrected Sherwood Number (Ref: Abramzon and Sirignano [1989])
        const scalar ShFilmCorr = 2 + (ShFilm - 2) / FM;

        //- Flux calculation    
        mFlux_ =  -(rhoFilm_ * mDiffFilm * ShFilmCorr)  / this->d() * log(1.0 + BM()); 

        //Heat Transfer due to radiation
        qR_ = cloud.fpRadiationModel().qR(TD_, mFlux_);  
         
        //- Internal heating
        qD_ = ((haGas_ - haSurf_) / (pow((1.0 + BM_), ZFilm_) - 1.0)) - (haSurf_ - hD_) + qR_;
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const label celli = this->cell();

    ParcelType::calc(cloud, td, dt); //-calc of KinematicParcel being called

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Gas phase mixture fraction calculation
    tetIndices tetIs = this->currentTetIndices();  
        
    //- Calculate mean and variance of mixture fraction.    
    scalar fRef = td.fInterp().interpolate(this->coordinates(), tetIs);  

    //- Bring in "cloudProperties" dictionary
    const dictionary popeSubModels(cloud.pope().cloudProperties().subDict("subModels"));
    const dictionary MMCcurlCoeffs(popeSubModels.subDict("MMCcurlCoeffs"));
    const dictionary Ximi(MMCcurlCoeffs.subDict("Ximi"));

    const scalar ri = readScalar(MMCcurlCoeffs.lookup("ri"));
    const scalar fm = readScalar(Ximi.lookup("fm"));


    //* Set hash tables for Multiphase - MMC coupling
    const labelHashTable& XiRIndexes = cloud.pope().mixing().XiR().rVarInXiR();
    const labelHashTable& XiCIndexes = cloud.pope().coupling().XiC().cVarInXiC();      

    YGas_.setSize(this->nThermoSpecies(), 0.0);
    YFilm_.setSize(this->nThermoSpecies(), 0.0);
    YSurf_.setSize(this->nThermoSpecies(), 0.0);
    YFuelVap_.setSize(cloud.nFuelSpecies(), 0.0);

    YFuelVap_ = YFuel_; // start with the same composition for Fuel Vapor and liquid fuel (does only change if diff evap is considered)

    // Initialize super cell stochastic particles              
    superCelli_ = cloud.pope().pManager().getSuperCellID(celli);    
    typedef typename basicAerosolReactingPopeCloud::particleType particleType;
    
    // ToDo: This operation can be expensive. Better to precompute the 
    //       particles in each super cell with 
    //       getParticlesInSuperCellList(cloud.pope()) and then hand over 
    //       the particles for this super cell
    DynamicList<particleType*> scParticles = 
        cloud.pope().pManager().getParticlesInSuperCell(cloud.pope(),superCelli_);
        
    
    std::random_shuffle(scParticles.begin(),scParticles.end());

    if(scParticles.size() <= 0.0)   
    {
        return;  
    }
    label id = 99999;  // far larger than the max number of sp in one super cell.
    scalar minDist = GREAT;
   
    vector posFP = this->position();

    forAll(scParticles, np)
    {
        vector posSP = scParticles[np]->position();

        vector rDist = posSP - posFP;
        scalar fDist = scParticles[np]->XiR()[XiRIndexes["f"]] - fRef;
        scalar dist = pow(rDist.x()/ri,2) + pow(rDist.y()/ri,2) + pow(rDist.z()/ri,2) + pow(fDist/fm,2);
        if(dist <= minDist)
        {
            minDist = dist; 
            id = np;
        }
    }
    
    fGas_ = scParticles[id]->XiC()[XiCIndexes["z"]];
    TGas_ = scParticles[id]->T();
    YGas_ = scParticles[id]->Y();
    haGas_ = scParticles[id]->hA();                             

    // initialize envelope model
    cloud.envelope().initialize(YGas_,TGas_,superCelli_);

    //- Access to surface mixtrue fraction
    this->particleH(cloud,td);
    calcSurfaceProperties(cloud,td);

    //- Access calc therompysical values
    calcThermophysicalValues(cloud,td); 

    //- Initial paramters for Crank Nicolson integration term f0  
    scalar np0 = this->nParticle();  
    scalar M0 = -mFlux_*this->areaS();        
    this->Cp_ = cloud.cpModel().Cp(this->TD_);
    scalar F0 = M0*qD_/( this->mass() * this->Cp() ); 
    scalar hDvol = haFuelVap_;

    // Sources
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    // Momentum transfer from the particle to the carrier phase for mass 
    vector dMTrans = vector::zero;

    // Mass rate from the particle to the carrier phase 
    scalar dMTransRate = 0.0; 
          
    scalar TSurfPrevious = TD_;

    // Reset surface state properties
    scalar dm_fp = 0.0;
    scalar timeFactor = 1.0; 

    scalar hGTD = cloud.pope().composition().particleMixture(scParticles[id]->Y()).Ha(scParticles[id]->pc(), TSurfPrevious);
    scalar hGTG = scParticles[id]->hA();

    //- Update new particle temperature and size~~~~~~~~~~~~~~~~~~~~~~~
    calcParticleUpdate(cloud,td, M0, F0, scParticles[id]->m(), hGTG, hGTD, hDvol, dm_fp, timeFactor, dt);

    //- update source terms for LES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    calcSourceLES(cloud, td, dt, dMTrans, dMTransRate, celli); 

    //- LES-LFP Coupling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //- Accumulate carrier phase source terms
    //- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    if (cloud.solution().coupled() && !noMassTransfer_)
    {           
        //- Update momentum transfer due to mass transfer 
        cloud.MTrans()[celli] +=  np0*dMTrans;  

        //- Update mass transfer rate for mixture fraction
        cloud.MTransRate()[celli] += np0*dMTransRate;
    }

    scalar dm_sp = np0 * dm_fp;   

    if(dm_sp > 0.0)
    {
        if(scParticles.size() > 0)
        {
            //- LFP --> MMC coupling ~~~~~~~~~~~~~~~~~~~~~~~~~~                   
            scalar m_sp = scParticles[id]->m();
            scalar weighting = dm_sp / (m_sp + dm_sp);
            scalar tauPhi = dt / (weighting + ROOTVSMALL);

            if(!noMassTransfer_) //- if fSurf < fGas no mass is transfered, but enthalpy still changes due to heating
            {                    //  (not good: the treatment of fSurf < fGas needs to be improved)
	            //.. Mixture fraction
                scalar c_z = 1.0; 
                scalar phi_z = scParticles[id]->XiC()[XiCIndexes["z"]];
                scParticles[id]->dpXiCsource()[XiCIndexes["z"]] = this->calcMMCSource(dt, c_z, phi_z, tauPhi);

	            //.. Species mass fraction
                forAll(cloud.linkPG(), nfs)
                {
                    scalar c_y = cloud.YProducts()[nfs];
                    scalar phi_y = scParticles[id]->Y()[cloud.linkPG()[nfs]];
                    scParticles[id]->dpYsource()[cloud.linkPG()[nfs]] = this->calcMMCSource(dt, c_y, phi_y, tauPhi);
                } 

                //.. Mass 
                scParticles[id]->dpMsource() = dm_sp; 
            }

            //...Enthalpy   
            scalar c_h = hD_ - timeFactor * (qD_ - qR_) + cloud.constProps().FH();
            scalar phi_h = scParticles[id]->hA(); 
            scParticles[id]->dp_hAsource() = this->calcMMCSource(dt, c_h, phi_h, tauPhi);

	        //.. Update stochastic particle properties                                     
            scParticles[id]->updateDispersedSources(); 
            scParticles[id]->wt()=scParticles[id]->m()/cloud.pope().deltaM();               

            scParticles[id]->T() = cloud.pope().composition().particleMixture(scParticles[id]->Y()).THa(scParticles[id]->hA(), scParticles[id]->pc(), scParticles[id]->T());
	    scParticles[id]->hEqv() = cloud.pope().composition().particleMixture(scParticles[id]->Y()).Hs(scParticles[id]->pc(),scParticles[id]->T());
        }
    }

    this->particleH(cloud,td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam:: sprayThermoParcel<ParcelType>::calcParticleUpdate
(
    TrackCloudType& cloud,
    trackingData& td,  
    scalar M0, 
    scalar F0, 
    scalar msp, 
    scalar& hGTG, 
    scalar& hGTD, 
    scalar& hDvol,
    scalar& dm_fp,
    scalar& timeFactor,
    const scalar dt
)
{
    scalar dmdt;
    scalar M1 = 0.0;
    scalar epsilon = (fSurf_ - fGas_) / max(fGas_, ROOTVSMALL);
    scalar threshold = 1e-7;
    if(epsilon < threshold)
    {
        mFlux_ = 0.0;
        dmdt = 0.0;
    }
    else 
    {
        scalar ratio = 0.0;
        M1 = -mFlux_ * this->areaS(this->d());
        dmdt = M0 * ratio + M1 * (1.0 - ratio);
    }

    //- update the accumulated volatile
    scalar cutoff = 1.0;
    scalar massPre = this->mass();
    if(dmdt == 0.0)
    {
        dm_fp = 0.0;
    }
    else if(dmdt < 0.0)
    {
        Pout<< " For fuel particle: " << this->origId() << " is having a negative dmdt (for volatile yield). " << nl;
        Pout<< " M0   = " << M0 << nl;        
        Pout<< " M1   = " << M1 << nl;
        Pout<< " dmdt = " << dmdt << nl;
        FatalErrorIn
        (
            " Wrong dmdt on FP "
        )   << exit(FatalError);           
    }
    else
    {
        //-Integrate for volatile mass 
        scalar mV = 0;
        scalar beta =0;
        const scalar deltaM = 
            cloud.MIntegrator().delta(mV, dt, dmdt, beta);
        mV = mV + deltaM;
        scalar massPre = this->mass();   
     
        if ( massPre > mV && !noMassTransfer_ )
        {
            //- update particle diameter and mass
            this->d() = pow((6.0 * (massPre - mV)/(pi * this->rho())), 1.0/3.0); 

	        //- update YFuel: change YFuel according to YFuelVap and transfered mass
	        forAll(YFuelVap_, nfs) YFuel_[nfs] = (YFuel_[nfs]*massPre - YFuelVap_[nfs]*mV) / (massPre - mV);
        }

        scalar massNew = this->mass();
        dm_fp = dmdt * dt; //massPre - massNew;

        if(dm_fp < 0.0)
        {
            Pout<< " M0   = " << M0 << nl;        
            Pout<< " M1   = " << M1 << nl;
            Pout<< " dmdt = " << dmdt << nl;
            Pout<< " mPre = " << massPre << nl;
            Pout<< " mNew = " << massNew << nl;
            Pout<< " dm   = " << dm_fp << nl;
            Pout<< " ratio= " << massPre / massNew << nl;
            
            FatalErrorIn
            (
                " Negative dm_fp on FP "
            )   << exit(FatalError); 
        }

    }    

    scalar F1 = 0.0;
    scalar hDTG = 0.0;
    scalar hDTD = 0.0;
    scalar dtMaxFp = GREAT;
    scalar dtMaxSp1 = GREAT;

    if(dm_fp > 0.0)
    {
        //- integrate to update TD_ 
        if(epsilon < threshold)
        {
            qD_ = 0.0;
            F0 = 0.0; // overwrite F0 to cut off the fuel particle heat transfer
        }   // Preventing FP exception in qD calculation
        else
        {
            this->Cp_ = cloud.cpModel().Cp(TD_);
            F1 = M1 * qD_ / (this->mass() * this->Cp());
        }
        
        scalar dm_sp = this->nParticle_ * dm_fp;
        scalar tauPhi = msp * dt / (dm_sp + ROOTVSMALL);
        scalar tauFp = massPre * dt / (dm_fp + ROOTVSMALL);

        scalar constantSp1 = fabs((hGTG - hGTD) / (hDvol - qD_ + qR_ - hGTG + ROOTVSMALL)); 

        if(constantSp1 < 1.0)
        {
            dtMaxSp1 = -tauPhi * log(1.0 - constantSp1);
        }
        
	    hDTG = cloud.pope().composition().particleMixture(this->YLiqToYGas(cloud,td,YFuelVap_)).Ha(pc_, TGas_) - cloud.constProps().LH();
        hDTD = hD_; // ok for now, but for diff evap hDTG = liq(YFuelVap@TD), not liq(YFuel@TD)?
        scalar constantFp = (hDTG - hDTD) / (qD_ + ROOTVSMALL);
        dtMaxFp = tauFp * constantFp;
        if(dtMaxFp < 0.0)
        {
            dtMaxFp = dt;
            // This deals with two situations: 
            // 1. When radiation is dominating in heat transfer. 
            // 2. When qD is ridiculously large due to a small BM. 
        }

        timeFactor = min(min(dtMaxSp1 / dt, dtMaxFp / dt), 1.0);

        if(dtMaxSp1 < dt || dtMaxFp < dt)
        {
            cloud.parcelIdAddOne(this->origId());
        }
    }

    scalar dTdt = timeFactor * cutoff * this->temperatureUpdateSwitch(F0, F1); // Preventing mass condensation
    scalar T_Dpre = TD_;  
    scalar beta =0.0;  
    const scalar deltaTD =
        cloud.TIntegrator().delta(TD_, dt, dTdt, beta);
    TD_ = min
          (
              max
              (
                  TD_ + deltaTD,
                  cloud.constProps().TMin()
              ),
              cloud.constProps().TMax()
          );     

    if(TD_ <= 200.0 || TD_ >= 3000.00)
    {
        Pout<< " -->Warning! Fuel particle: " << this->origId() << " wrong temperature " <<nl; 
        Pout<< " TG   = " << TGas_ << nl;
        Pout<< " TD   = " << T_Dpre << nl;
        Pout<< " hDTG = " << hDTG << nl;
        Pout<< " hDTD = " << hDTD << nl;
        Pout<< " qD   = " << qD_ << nl;       
        Pout<< " dm_fp     = " << dm_fp << nl;
        Pout<< " mD        = " << massPre << nl;
        Pout<< " dt*fp     = " << dtMaxFp << nl;
        Pout<< " dt        = " << dt << nl;
        Pout<< " timeFactor = " << timeFactor << nl;
        FatalErrorIn
        (
            " Wrong temperature on FP "
        )   << exit(FatalError);   
    }

    //- update hD
    this->particleH(cloud,td);

    //- this member function will also return the cutoff ratio for backward couplings
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::sprayThermoParcel<ParcelType>::calcSourceLES
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    vector& dMTrans,
    scalar& dMTransRate,                
    const label celli
)
{
    //-Update momentum source term for mass 
    dMTrans += - dt* mFlux_*this->areaS(this->d())* this->U(); 

    //-Update continuity and mixture fraction source trerm for mass flow rate
    dMTransRate += - dt* mFlux_*this->areaS(this->d()); 
}


//- turns composition of n_fuel components into n_thermo components 
template<class ParcelType>
template<class TrackCloudType>
Foam::scalarField Foam::sprayThermoParcel<ParcelType>::YLiqToYGas
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalarField YLiq 
)
{
    scalarField YGas;
    YGas.setSize(this->nThermoSpecies(), 0.0);
    scalar sumY = 0.0;

    forAll(cloud.linkFG(), nfs)
    {
        YGas[cloud.linkFG()[nfs]] = YLiq[nfs];
    }

    forAll(YGas,i)
    {
        sumY += YGas[i];
    }
    YGas /= sumY;
    return YGas;
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::sprayThermoParcel<ParticleType>::getStatisticalData
(
    const wordList& vars
) const
{
    // Create with 6 scalar fields for the mandatory entries
    DynamicList<scalar> container(6);

    container.append(this->mesh().time().time().value());
    container.append(this->position().x());
    container.append(this->position().y());
    container.append(this->position().z());

    // if no word list is given all variables are sampled
    if (vars.empty())
    {
        table_.storeAllVars(container);
    }
    else
    {
        table_.storeVarsByList(container,vars);
    }

    return container;
}


template<class ParticleType>
Foam::DynamicList<Foam::word>
Foam::sprayThermoParcel<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    DynamicList<word> varNames(6);
    
    varNames.append("time");
    varNames.append("px");
    varNames.append("py");
    varNames.append("pz");


    // if no word list is given all variables are sampled
    if (vars.empty())
        varNames.append(table_.getAllVarNames());
    else
        varNames.append(vars);
    
    return varNames;
}


template<class ParticleType>
void Foam::sprayThermoParcel<ParticleType>::initStatisticalSampling()
{
    // sprayThermoParcel does not derive from ItoPopeParticle and 
    // therefore we need to add here the position and time

    this->nameVariableLookUpTable().addNamedVariable
        ("D",this->d_);

    this->nameVariableLookUpTable().addNamedVariable
        ("ux",this->U_.x());
    this->nameVariableLookUpTable().addNamedVariable
        ("uy",this->U_.y());
    this->nameVariableLookUpTable().addNamedVariable
        ("uz",this->U_.z());
    
    this->nameVariableLookUpTable().addNamedVariable
        ("nParticle",this->nParticle_);
    this->nameVariableLookUpTable().addNamedVariable
        ("fGas",fGas_);
    this->nameVariableLookUpTable().addNamedVariable
        ("fSurf",fSurf_);
    this->nameVariableLookUpTable().addNamedVariable
        ("TD",TD_);
    this->nameVariableLookUpTable().addNamedVariable
        ("TGas",TGas_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const sprayThermoParcel<ParcelType>& p
)
:
    ParcelType(p),
    TD_(p.TD_),
    TGas_(p.TGas_),
    TFilm_(p.TFilm_),
    TSurf_(p.TSurf_),
    pc_(p.pc_),
    fGas_(p.fGas_),
    fSurf_(p.fSurf_),
    fFilm_(p.fFilm_),
    nThermoSpecies_(p.nThermoSpecies_),
    nFuelSpecies_(p.nFuelSpecies_),
    YFuel_(p.YFuel_),
    YFuelVap_(p.YFuelVap_),
    YGas_(p.YGas_),
    YFilm_(p.YFilm_),
    YSurf_(p.YSurf_),
    haGas_(p.haGas_),
    hD_(p.hD_),
    haSurf_(p.haSurf_),
    haFilm_(p.haFilm_),
    haFuelVap_(p.haFuelVap_),
    muFilm_(p.muFilm_),
    rhoFilm_(p.rhoFilm_),
    qR_(p.qR_),
    qD_(p.qD_),
    mFlux_(p.mFlux_),
    molWtSstate_(p.molWtSstate_),   
    ZFilm_(p.ZFilm_),
    BM_(p.BM_)
{
    initStatisticalSampling();
}


template<class ParcelType>
Foam::sprayThermoParcel<ParcelType>::sprayThermoParcel
(
    const sprayThermoParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    TD_(p.TD_),
    TGas_(p.TGas_),
    TFilm_(p.TFilm_),
    TSurf_(p.TSurf_),
    pc_(p.pc_),
    fGas_(p.fGas_),
    fSurf_(p.fSurf_),
    fFilm_(p.fFilm_),
    nThermoSpecies_(p.nThermoSpecies_),
    nFuelSpecies_(p.nFuelSpecies_),
    YFuel_(p.YFuel_),
    YFuelVap_(p.YFuelVap_),
    YGas_(p.YGas_),
    YFilm_(p.YFilm_),
    YSurf_(p.YSurf_),
    haGas_(p.haGas_),
    hD_(p.hD_),
    haSurf_(p.haSurf_),
    haFilm_(p.haFilm_),
    haFuelVap_(p.haFuelVap_),
    muFilm_(p.muFilm_),
    rhoFilm_(p.rhoFilm_),
    qR_(p.qR_),
    qD_(p.qD_),
    mFlux_(p.mFlux_),
    molWtSstate_(p.molWtSstate_),   
    ZFilm_(p.ZFilm_),
    BM_(p.BM_)
{
    initStatisticalSampling();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "sprayThermoParcelIO.C"

// ************************************************************************* //
