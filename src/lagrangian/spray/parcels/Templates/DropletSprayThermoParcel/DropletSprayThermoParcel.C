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
#include "physicoChemicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    tetIndices tetIs = this->currentTetIndices();

    td.Cpc() = td.CpInterp().interpolate(this->coordinates(), tetIs);

    td.Tc() = td.TInterp().interpolate(this->coordinates(), tetIs);

    if (td.Tc() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed temperature in cell " << this->cell()
                << " to " << cloud.constProps().TMin() <<  nl << endl;
        }

        td.Tc() = cloud.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const label celli = this->cell();
    const scalar massCell = this->massCell(td);

    td.Uc() += cloud.UTrans()[celli]/massCell;

//    tetIndices tetIs = this->currentTetIndices();
//    Tc_ = td.TInterp().interpolate(this->coordinates(), tetIs);


    // ==> IMPORTANT <==
    // Not sure if this makes sense here, as all thermo properties
    // are coupled back through the mmcCloud
    //
    // const scalar CpMean = td.CpInterp().psi()[celli];
    // td.Tc() += cloud.hsTrans()[celli]/(CpMean*massCell);

    // if (td.Tc() < cloud.constProps().TMin())
    // {
    //     if (debug)
    //     {
    //         WarningInFunction
    //             << "Limiting observed temperature in cell " << celli
    //             << " to " << cloud.constProps().TMin() <<  nl << endl;
    //     }

    //     td.Tc() = cloud.constProps().TMin();
    // }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::calcReferenceValues
(
    TrackCloudType& cloud,
    const scalar pInf,          // Infinity pressure
    const scalar Ts,            // Surface temperature
    const scalar TInf,          // Infinity temperature
    const List<scalar>& Ys,     // Surface species composition
    const List<scalar>& YInf,   // Infinity species composition
    const label jFuel,       // Index of the fuel mass fraction in Y
    scalar& TRef,               // Reference temperature
    scalar& rhoRef,             // Density at reference conditions
    scalar& muRef,              // Viscosity at reference conditions
    scalar& kappaRef,           // heat conductivity at reference condition
    scalar& CpRef,              // heat capacity at reference conditions
    scalar& CpFRef          // heat capacity at reference conditions (pure liquid mass fraction)
) const
{
    // The pressure value is selected from the 

    // Calculate the reference temperature and mass fraction
    TRef = (TInf + 2.0*Ts)/3.0;
    List<scalar> YRef(Ys);
    
    forAll(YRef,i)
    {
        YRef[i] = (YInf[i] + 2.0*Ys[i])/3.0;
    }

    // Determine the thermophysical properties 
    auto& thermoModels = cloud.thermoModels();    

    // Update with the current species mass fractions to get the correct
    // mixture properties
    thermoModels.update(YRef);

    rhoRef = thermoModels.rho(pInf,TRef);
    muRef = thermoModels.mu(pInf,TRef);
    kappaRef = thermoModels.kappa(pInf,TRef);
    CpRef = thermoModels.Cp(pInf,TRef);

    CpFRef = thermoModels.Cp(pInf,TRef,jFuel);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Do not calculate for very small time steps, otherwise numerical error
    if (dt < 1E-15)
        return;

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar mD0 = this->mass0_;
    const scalar T0 = this->T_;

    // Infinity conditions
    scalar TInf=0;
    scalar pInf=0;
    scalarList YInf;

    // Get the stochastic particles to determine the infinity conditions
    auto& couplingModel = cloud.dropletToMMCCoupling();
    couplingModel.getGasPhaseConditions
    (
        *this,
        td,
        cloud.mmcCloud(),
        TInf,
        pInf,
        YInf
    );

    // If the species index is unset, set it now 
    if (jFuel_ == -1)
    {
        const word& speciesName = cloud.dropletSpeciesName();
        jFuel_ = cloud.thermo().carrier().species()[speciesName];
    }

    //- Dummy values that are set in the time integration
    scalar muRef = SMALL;
    scalar rhoRef = SMALL;
    scalar mDNew = SMALL;
    scalar TDNew = SMALL;


    switch (cloud.integrationScheme())
    {
        case DropletSprayThermo::integrationScheme::EulerExplicit:
            timeIntegrationEulerExplicit
            (
                cloud,
                td,
                muRef,
                rhoRef,
                mDNew,
                TDNew,
                pInf,
                this->T_,
                TInf,
                YInf,
                dt
            );
            break;

        case DropletSprayThermo::integrationScheme::CrankNicolson:
            timeIntegrationCrankNicolson
            (
                cloud,
                td,
                muRef,
                rhoRef,
                mDNew,
                TDNew,
                pInf,
                this->T_,
                TInf,
                YInf,
                dt
            );
            break;

        default:
            FatalError 
                << "No valid integration scheme is specified for droplet evaporation"
                << exit(FatalError);
            break;
    }

    // Update droplet diameter and density
    this->T_ = TDNew;
    this->mass0_ = mDNew;
    this->rho() = cloud.liquidProperties().rhoL(pInf,this->T_);
    this->d_ = std::pow(mDNew*6.0/(this->rho()*pi),1.0/3.0);
    this->Cp_ = cloud.liquidProperties().CpL(pInf,this->T_); //Changed 
    

    if (td.keepParticle == false)
        this->active_ = false;

    if (cloud.solution().coupled() && this->active())
    {
        // Couple the source terms to the mmcCloud
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Calculate the total mass and heat transfer
        const scalar totalMDot = -np0*(mDNew-mD0)/dt;

        // Get total time step
        const scalar trackTime = cloud.solution().trackTime();

        mDot_+=-np0*(mDNew-mD0)/trackTime;

        // Calculate total enthalpy exchange with the gas
        const scalar hA_new = cloud.liquidProperties().Ha(pInf,this->T_);
        const scalar hA_old = cloud.liquidProperties().Ha(pInf,T0);

        const scalar hDot = -np0*(mDNew*hA_new - mD0*hA_old)/dt;
        couplingModel.coupleSources
        (
            *this,
            td,
            cloud.mmcCloud(),
            totalMDot,
            hDot,
            jFuel_,
            dt
        );

        // Couple the source terms to the Eulerian field
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Explicit momentum source for particle
        vector Su = Zero;
    
        // Linearised momentum source coefficient
        scalar Spu = 0.0;
    
        // Momentum transfer from the particle to the carrier phase
        vector dUTrans = Zero;

        // Calcualte the Reynoldsnumber
        const scalar Re = this->Re(rhoRef, this->U_, td.Uc(), this->d_, muRef);

        // Store old particle velocity
        const vector uD0 = this->U_;

        // Calculate new particle velocity
        this->U_ =
            this->calcVelocity(cloud, td, dt, Re, muRef, mD0, Su, dUTrans, Spu);


        //  Accumulate carrier phase source terms
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;

        // Add the evaporation term:
        cloud.UTrans()[this->cell()] += -np0*uD0*(mDNew-mD0);

        // Update momentum transfer coefficient
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update mass transfer term
        cloud.rhoTrans()[this->cell()] += totalMDot*dt;
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::timeIntegrationEulerExplicit
(
    TrackCloudType& cloud,
    trackingData& td,
    scalar& muRef,              // Viscosity at reference conditions
    scalar& rhoRef,             // Gas density at reference conditions
    scalar& mDNew,              // Mass of the droplet after evaporation
    scalar& TDNew,              // Droplet temperature after evaporation
    const scalar pInf,          // Infinity pressure
    const scalar Ts,            // Surface temperature
    const scalar TInf,          // Infinity temperature
    const List<scalar>& YInf,   // Infinity species composition
    const scalar dt             // Time step
) const
{
    // Initial values at start point
    const scalar& TD0 = this->T_;
    const scalar mD0 = this->mass0_;
    const scalar TSat = cloud.liquidProperties().TSat(pInf)-1E-3;
    const scalar CpL = cloud.liquidProperties().CpL(pInf,this->T_);
    scalar mDot, QDot;


    // Calculate the mass and heat exchange
    mDotAndQDot(cloud,td,mDot,QDot,muRef,rhoRef,this->d_,pInf,this->T_,TInf,YInf,jFuel_);

    mDNew = mD0;

    // Change in mass
    const scalar deltaM = mDot*dt;
    if (deltaM < mD0)
    {
        // Get Cp value of the fuel at pInf and T droplet
        mDNew = mD0-deltaM;
        TDNew = std::min(TD0 + (QDot*dt/(CpL*mD0)),TSat);
    }
    else
    {
        mDNew = 0;
        td.keepParticle = false;
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::timeIntegrationCrankNicolson
(
    TrackCloudType& cloud,
    trackingData& td,
    scalar& muRef,              // Viscosity at reference conditions
    scalar& rhoRef,             // Gas density at reference conditions
    scalar& mDNew,              // Mass of the droplet after evaporation
    scalar& TDNew,              // Droplet temperature after evaporation
    const scalar pInf,          // Infinity pressure
    const scalar Ts,            // Surface temperature
    const scalar TInf,          // Infinity temperature
    const List<scalar>& YInf,   // Infinity species composition
    const scalar dt             // Time step
) const
{
    // Define some parameters controlling the convergence
    const unsigned int iterMax = 50; // Maximum number of NR iterations
    const scalar eps   = 1e-12;  // Maximum residual for NR iteration (delta)
    const scalar eps_m = 1e-03;  // Relative difference of m for calculating Jacobian
    const scalar eps_T = 1e-03;  // Relative difference of T for calculating Jacobian

    const scalar TSat = cloud.liquidProperties().TSat(pInf)-1E-3;
    const scalar CpL0 = cloud.liquidProperties().CpL(pInf,this->T_);
    scalar rhoL = this->rho();


    // Initial values at start point
    const scalar& TD0 = this->T_;
    const scalar& dD0 = this->d_;
    const scalar mD0 = this->mass0_;

    // Give default value to new state
    mDNew = mD0;
    TDNew = TD0;
    scalar dDNew = dD0;

    // Dummy values
    scalar mDot, QDot;

    // Calculate the mass and heat exchange based on current data
    mDotAndQDot(cloud,td,mDot,QDot,muRef,rhoRef,dD0,pInf,TD0,TInf,YInf,jFuel_);
    const scalar rhs_m_o  = -mDot;
    const scalar rhs_T_o  =  QDot/(CpL0*mDNew);

    // Estimated evaporation time
    scalar tauEvap = 1000.0*dt;
    if (mDot > 0.0)
    {
        tauEvap = 1.5*mD0/(mDot);
    }

    if (tauEvap < 2.0*dt)
    {
        mDNew = 0;
        td.keepParticle = false;
        return;
    }

    const scalar mRef = mD0;  // Reference values for normalization
    const scalar TRef = TD0;
    unsigned int iter = 0;

    // Update the heat capacity based on the current new temperature
    scalar CpL = cloud.liquidProperties().CpL(pInf,TDNew);

    for (iter=0; iter < iterMax; ++iter)
    {
        CpL = cloud.liquidProperties().CpL(pInf,TDNew);
        // get the heat and mass transfer rates
        mDotAndQDot(cloud,td,mDot,QDot,muRef,rhoRef,dDNew,pInf,TDNew,TInf,YInf,jFuel_);
        const scalar rhs_m_n  = -mDot;
        const scalar rhs_T_n  = QDot/(CpL*mDNew);

        const scalar Fm = (mDNew-mD0-0.5*dt*(rhs_m_o+rhs_m_n))/mRef;
        const scalar FT = (TDNew-TD0-0.5*dt*(rhs_T_o+rhs_T_n))/TRef;
        
        // Calculate components of the Jacobian matrix
        const scalar dDEps = pow(6.0/pi*(mDNew+(eps_m*mRef))/rhoL,1.0/3.0);
        mDotAndQDot(cloud,td,mDot,QDot,muRef,rhoRef,dDEps,pInf,TDNew,TInf,YInf,jFuel_);
        const scalar rhs_m_mp  = -mDot;
        const scalar rhs_T_mp  = QDot/(CpL*mDNew);

        const scalar TEps = TDNew-(eps_T*TRef);
        const scalar CpL_TEps = cloud.liquidProperties().CpL(pInf,TEps);
        const scalar rhoL_TEps = cloud.liquidProperties().rhoL(pInf,TDNew);
        const scalar dD_TEps = std::pow(mDNew*6.0/(rhoL_TEps*pi),1.0/3.0);
        mDotAndQDot(cloud,td,mDot,QDot,muRef,rhoRef,dD_TEps,pInf,TEps,TInf,YInf,jFuel_);
        const scalar rhs_m_Tm  = -mDot;
        const scalar rhs_T_Tm  = QDot/(CpL_TEps*mDNew);

        const scalar dFmdm = 1.0-(0.5*dt*(rhs_m_mp-rhs_m_n)/mRef)/eps_m;
        const scalar dFmdT = 0.0-(0.5*dt*(rhs_m_n-rhs_m_Tm)/mRef)/eps_T;
        const scalar dFTdm = 0.0-(0.5*dt*(rhs_T_mp-rhs_T_n)/TRef)/eps_m;
        const scalar dFTdT = 1.0-(0.5*dt*(rhs_T_n-rhs_T_Tm)/TRef)/eps_T;
        
        // Solve for new state
        const Pair<scalar> delta = matrixSolver2By2(dFmdm,dFmdT,dFTdm,dFTdT,-Fm,-FT);

        // ~ Final Step ~

        // Avoid heating beyond boiling point
        mDNew = mDNew+delta[0]*mRef;
        TDNew = std::min(TDNew+delta[1]*TRef,TSat);
        rhoL = cloud.liquidProperties().rhoL(pInf,TDNew);
        dDNew = std::pow(mDNew*6.0/(rhoL*pi),1.0/3.0);

        // Check if calculated values are physical
        if (mDNew <=0)
        {
            mDNew = 0;
            td.keepParticle = false;
            break;
        }
        
        // Residual control
        const scalar res = sqrt((pow(delta[0],2.0)+pow(delta[1],2.0))/2.0); // L2 norm
        if (res < eps) 
            break;
        
    }
    if (iter==iterMax)
    {
        Pout << "Warning: Crank-Nicolson method did not converge within " << iter << " iterations!\n";
        // Use explicit Euler method to evolve droplet properties in time
        mDNew = mD0+dt*rhs_m_o;
        TDNew = TD0+dt*rhs_T_o;
    }
}


template<class ParcelType>
Foam::Pair<Foam::scalar>
Foam::DropletSprayThermoParcel<ParcelType>::matrixSolver2By2
(
    const scalar a11,
    const scalar a12,
    const scalar a21,
    const scalar a22,
    const scalar b1,
    const scalar b2
) const
{
    // Solve a system of linear equations of the form A*x=b with dimension two
    const scalar detA = a11*a22-a12*a21;
    Pair<scalar> x;

    x[0] = -(a12*b2-a22*b1)/detA;
    x[1] =  (a11*b2-a21*b1)/detA;
    return x;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DropletSprayThermoParcel<ParcelType>::mDotAndQDot
(
    TrackCloudType& cloud,
    trackingData& td,
    scalar& mDotF,              // Evaporated mass of fuel [kg/s]
    scalar& QDot,               // Heat flux [J/s]
    scalar& muRef,              // Viscosity at reference conditions
    scalar& rhoRef,             // Gas density at reference conditions
    const scalar dD,            // Droplet diameter
    const scalar pInf,          // Infinity pressure
    const scalar Ts,            // Surface temperature
    const scalar TInf,          // Infinity temperature
    const List<scalar>& YInf,   // Infinity species composition
    const label jFuel           // Index of the fuel mass fraction in Y
) const
{
    // Calculate the mole fraction of the fuel species at the droplet surface
    const scalar pSat = cloud.liquidProperties().pSat(Ts);

    // Mole fraction of the fuel at the surface
    const scalar XsFuel = min(pSat/pInf,1.0-1E-08);

    auto& thermoModels = cloud.thermoModels();
    // Get the molecular weight vector [g/mol]
    List<scalar> MW(YInf.size());
    forAll(MW,i)
    {
        MW[i] = thermoModels.MW(i);
    }

    // Calculate the mole fraction of all species at infinity conditions
    const List<scalar> XInf = YtoX(YInf,MW);

    // Check if XInf is fully saturated
    // If already saturated only heat the droplet
    if (XInf[jFuel] >= 1.0)
    {
        mDotF = 0;
        // Get the heat conductivity of the fuel
        const scalar kappa = cloud.thermo().carrier().kappa(pInf,Ts,jFuel);
        muRef = cloud.thermo().carrier().mu(pInf,Ts,jFuel);
        rhoRef = cloud.thermo().carrier().rho(pInf,Ts,jFuel);
        QDot = 2.0*pi*dD*kappa*(TInf-Ts);
        return;
    }

    // Calculate the mole fractions at the droplet surface, except fuel
    List<scalar> Xs(XInf.size());
    
    forAll(Xs,i)
        Xs[i] = XInf[i]*(1.0-XsFuel)/(1.0-XInf[jFuel]);
    Xs[jFuel] = XsFuel;

    // Now calculate again the mass fractions based on the mole fractions
    List<scalar> Ys = XtoY(Xs,MW);

    // Calculate the density at infinity conditions:
    thermoModels.update(YInf);    

    // Evaluate the average gas properties
    scalar TRef;    // reference temperature
    scalar kappaRef;    // reference heat conductivity
    scalar CpRef;   // heat capacity
    scalar CpFRef;  // heat capacity of the fuel species

    calcReferenceValues
    (
        cloud,
        pInf, Ts, TInf, Ys, YInf, jFuel,
        TRef, rhoRef, muRef, kappaRef, CpRef, CpFRef
    );

    // Mol weight of the gas
    scalar MW_Mix = 0;
    forAll(MW,i)
    {
        MW_Mix += MW[i]*XInf[i];
    }
    scalar Df = cloud.liquidProperties().Df(pInf,TRef,MW_Mix);

    // Calculate the non-dimensional numbers
    const scalar Sc = muRef/(rhoRef*Df);
    const scalar Pr = CpRef*muRef/(kappaRef);
    const scalar Le = Sc/Pr;
    const scalar Re = this->Re(rhoRef,this->U_,td.Uc(),dD,muRef);
    const scalar Sh = 1.0+pow(1.0+Re*Sc,1.0/3.0)*max(1.0,pow(Re,0.077));
    const scalar Nu = 1.0+pow(1.0+Re*Pr,1.0/3.0)*max(1.0,pow(Re,0.077));


    // Spalding mass transfer number
    const scalar BM  = (Ys[jFuel]-YInf[jFuel])/(1.0-Ys[jFuel]);
    
    const scalar FM  = pow(1.0+BM,0.7)*log(1.0+BM)/BM;
    const scalar ShMod = 2.0+(Sh-2.0)/FM;
    mDotF = pi*dD*ShMod*rhoRef*Df*log(1.0+BM);
    scalar BT = SMALL;
    if (fabs(BM)>1e-12)
    {
        scalar phi = (CpFRef/CpRef)*(Sh/Nu)*(1.0/Le);  // initial guess
        const scalar BT0 = pow(1.0+BM,phi)-1.0;
        BT = BT0;

        unsigned int iter = 0;
        const unsigned int iterMax = 100;
        const scalar eps = 1.0e-08;
        for (iter=0; iter<iterMax; ++iter)
        {
            const scalar FT = pow(1.0+BT,0.7)*log(1.0+BT)/BT;
            const scalar NuMod = 2.0+(Nu-2.0)/FT;
            phi = (CpFRef/CpRef)*(ShMod/NuMod)*(1.0/Le);
            const scalar BTNew = pow(1.0+BM,phi)-1.0;
            const scalar res = fabs(BTNew-BT);
            BT = BTNew;
            if (res < eps) 
                break;
        }
        if (iter == iterMax-1)
        {
            Pout << "Warning: Iteration of BT did not converge within " << iter << " iterations!\n";
            BT = BT0;
        }
        
        const scalar Lv = cloud.liquidProperties().Lv(pInf,Ts);
        QDot  = mDotF*(CpFRef*(TInf-Ts)/BT-Lv);
    }
    else
    {
        QDot = 2.0*pi*dD*kappaRef*(TInf-Ts);  // pure heat transfer
    }
}


template<class ParcelType>
Foam::List<Foam::scalar> 
Foam::DropletSprayThermoParcel<ParcelType>::YtoX
(
    const List<scalar>& Y,
    const List<scalar>& MW
) const
{
    // Calculate temporary mean weighted mole mass
    scalar meanYMW = 0;
    forAll(Y,i)
    {
        meanYMW += Y[i]/MW[i];
    }

    List<scalar> Xi(Y.size());

    forAll(Xi,i)
    {
        Xi[i] = (Y[i]/MW[i])/meanYMW;
    }

    return Xi;
}


template<class ParcelType>
Foam::List<Foam::scalar> 
Foam::DropletSprayThermoParcel<ParcelType>::XtoY
(
    const List<scalar>& X,
    const List<scalar>& MW
) const
{
    // Calculate temporary mean weighted mole mass
    scalar meanXMW = 0;
    forAll(X,i)
    {
        meanXMW += X[i]*MW[i];
    }

    List<scalar> Y(X.size());

    forAll(Y,i)
    {
        Y[i] = X[i]*MW[i]/meanXMW;
    }

    return Y;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DropletSprayThermoParcel<ParcelType>::DropletSprayThermoParcel
(
    const DropletSprayThermoParcel<ParcelType>& p
)
:
    ParcelType(p),
    T_(p.T_),
    Cp_(p.Cp_),
    jFuel_(p.jFuel_)
{
    initStatisticalSampling();
}


template<class ParcelType>
Foam::DropletSprayThermoParcel<ParcelType>::DropletSprayThermoParcel
(
    const DropletSprayThermoParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    T_(p.T_),
    Cp_(p.Cp_),
    jFuel_(p.jFuel_)
{
    initStatisticalSampling();
}

// * * * * * * * * * * * * Statistica Data Sampling * * * * * * * * * * * * * *

template<class ParticleType>
void Foam::DropletSprayThermoParcel<ParticleType>::initStatisticalSampling()
{
    nameVariableLookUpTable().addNamedVariable("mDot",mDot_);
    nameVariableLookUpTable().addNamedVariable("d",this->d_);
    nameVariableLookUpTable().addNamedVariable("T",T_);
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::DropletSprayThermoParcel<ParticleType>::getStatisticalData
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
Foam::DropletSprayThermoParcel<ParticleType>::getStatisticalDataNames
(
    const wordList& vars
) const
{
    DynamicList<word> varNames(6);
    
    varNames.append("time");
    // Use px, py, and pz instead of x,y,z to avoid confusen with the 
    // conditional variable z
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

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "DropletSprayThermoParcelIO.C"

// ************************************************************************* //
