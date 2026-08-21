/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * NGDE - Aerosol Dynamics solver * * * * * * * * * * * * //

// The following code is an implementation of the journal publication and supplementary
// material published in:
//
// A. Prakash, A.P. Bapat and M.R. Zachariah, "A Simple Numerical Algorithm and Software
// for Solution of Nucleation, Surface Growth, and Coagulation Problems", Aerosol Science 
// and Technology, 37:892-898, 2003. DOI: 10.1080/02786820390225826
// 
// This C code is used to solve for nucleation, coagulation and surface growth problems. 
// This example problem solves for characteristics of an Aluminum aerosol. However, it can 
// be used for any other material whose properties are known. 
//
// The solution algorithm involves a new approach to solve the Aerosol GDE, where the 
// particle volume space is divided into nodes of zero width. Using nodes allows us to 
// cover the entire size range of aerosol (1 nm to 10 micrometer) using just 40 nodes.
//
// Gregor Neuber (gregor.neuber@itv.uni-stuttgart.de), May 2017

// * * * * * * * * * * * * * * * * Variables * * * * * * * * * * * * * * * * * //

// v_[i]                  volume of particles at ith node (m^3/m^3 of aerosol)                   
// particle.psdProperties()    number of particles at ith node (#/m^3 of aerosol)                    
// X_                     splitting operator for coagulation (dimensionless)                   
// beta[i][j]             collision frequency function (m^3/s)                                   
// N1s[i]                 saturation monomer concentration over particle of size i (#/m^3)       
// coolrate_              cooling rate (K/s)                                                
// addterm,subterm        addition and subtraction terms in Nk due to surface growth (#/(m^3s))  
// sigma                  surface Tension (N/m)                                                
// Ps                     saturation Pressure (Pa)                                             
// ns                     monomer concentration at saturation (#/m^3 of aerosol)                 
// kstar                  critical node size (dimensionless)                                   
// dpstar                 critical particle size (m)                                           
// vstar                  volume of particles in the critical bin (m^3/m^3 of aerosol)           
// s1                     surface area of monomer (m^2)                                         
// S                      saturation ratio (dimensionless)                                     
// m1                     monomer mass (kg)                                                     
// zeta[i]                splitting operator for nucleation at the ith node (dimensionless)    

// * * * * * * * * * * * * * * * * Includes * * * * * * * * * * * * * * * * * //

#include "SectionalDropletSynthesis.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SectionalDropletSynthesis<CloudType>::SectionalDropletSynthesis
(
    const dictionary& dict,
    CloudType& owner
)
:
    SynthesisModel<CloudType>(dict,owner,typeName),
    
    speciesThermo_
    (
      dynamic_cast<const reactingMixture<ThermoType>&> (owner.thermo()).speciesData()
    ),
    
    saturation_(0),
    
    nucleationRate_(1),
    
    surfaceGrowthRate_(2),
    
    nPhysicalAerosolProperties_(3),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
    
    d0_(readScalar(this->coeffDict().lookup("d0"))),
       
    growth_(readScalar(this->coeffDict().lookup("growth"))),
            
    d_(0.0),
    
    a_(0.0),
    
    v_(0.0),
    
    X_(0.0),

    condensingSpeciesName_(this->coeffDict().lookup("condensingSpeciesName")),
    
    condensingSpeciesIndex_(0),
    
    liquidDensity_(readScalar(this->coeffDict().lookup("liquidDensity"))),
    
    coagulationModel_(this->coeffDict().lookup("coagulationModel")),
    
    coolrate_(readScalar(this->coeffDict().lookup("coolrate"))),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    kb_(Ru_/NA_)
{
    forAll(speciesThermo_, specieI)
    {
        if (owner.slgThermo().carrier().species()[specieI] == condensingSpeciesName_)
        {
            condensingSpeciesIndex_ = specieI;
        }
    }
        
    Info << "Condensing species: "        << owner.slgThermo().carrier().species()[condensingSpeciesIndex_]  << endl;
    Info << "Diameter of first section: " << d0_     << endl;
    Info << "Diameter growth: "           << growth_ << endl;
    Info << "Avogadro-Constant = "        << NA_     << endl;
    Info << "Gas constant = "             << Ru_     << endl;
    Info << "Boltzman constant = "        << kb_     << endl;
    
    Info << coagulationModel_ << " model is used as the coagulation model." << nl << endl;
        
    setSectionProperties();
    
    setSizeSplittingOperator();
}


template <class CloudType>
Foam::SectionalDropletSynthesis<CloudType>::SectionalDropletSynthesis
(
    const SectionalDropletSynthesis<CloudType>& cm
)
:
    SynthesisModel<CloudType>(cm),
    
    speciesThermo_(cm.speciesThermo_),
    
    saturation_(0),
    
    nucleationRate_(1),
    
    surfaceGrowthRate_(2),
    
    nPhysicalAerosolProperties_(3),
    
    nSections_(readLabel(this->coeffDict().lookup("nSections"))),
    
    d0_(readScalar(this->coeffDict().lookup("d0"))),
       
    growth_(readScalar(this->coeffDict().lookup("growth"))),
    
    d_(0.0),
    
    a_(0.0),
    
    v_(0.0),
    
    X_(0.0),
            
    condensingSpeciesName_(this->coeffDict().lookup("condensingSpeciesName")),
    
    condensingSpeciesIndex_(0),
    
    liquidDensity_(readScalar(this->coeffDict().lookup("liquidDensity"))),
    
    coagulationModel_(this->coeffDict().lookup("coagulationModel")),
    
    coolrate_(readScalar(this->coeffDict().lookup("coolrate"))),
    
    NA_(physicoChemical::NA.value()*1000),
    
    Ru_(physicoChemical::R.value()*1000),
    
    kb_(Ru_/NA_)
{
    forAll(speciesThermo_, specieI)
    {
        if (cm.slgThermo().carrier().species()[specieI] == condensingSpeciesName_)
        {
            condensingSpeciesIndex_ = specieI;
        }
    }

    Info << "Condensing species: "        << cm.slgThermo().carrier().species()[condensingSpeciesIndex_]  << endl;
    Info << "Diameter of first section: " << d0_     << endl;
    Info << "Diameter growth: "           << growth_ << endl;
    Info << "Avogadro-Constant = "        << NA_     << endl;
    Info << "Gas constant = "             << Ru_     << endl;
    Info << "Boltzman constant = "        << kb_     << endl;
    
    Info << coagulationModel_ << " model is used as the coagulation model." << nl << endl;
        
    setSectionProperties();
    
    setSizeSplittingOperator();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SectionalDropletSynthesis<CloudType>::~SectionalDropletSynthesis()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::synthesize()
{
    Info << "Solid particle inception..." << endl;    
    
    forAllIters(this->owner(), iter)
    {
        particleInception(iter());
    }
    
    Info << "Solid particle surface growth..." << endl;    
    
    forAllIters(this->owner(), iter)
    {
        particleSurfaceGrowth(iter());
    }
    
    Info << "Solid particle coagulation..." << endl;    
    
    forAllIters(this->owner(), iter)
    {
        particleCoagulation(iter());
    }
    
    Info << "Gas phase cooling..." << endl;    
    
    forAllIters(this->owner(), iter)
    {
        gasPhaseCooling(iter());
    }
    
    Info << "Finished!" << endl;    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::setSectionProperties()
{
    d_.setSize(nSections_+1, 0.0);
    a_.setSize(nSections_+1, 0.0);
    v_.setSize(nSections_+1, 0.0);

    v_[0] = speciesThermo_[condensingSpeciesIndex_].W()/(liquidDensity_*NA_);
    
    // calculation of geometric spacing factor which depends on the number of nodes
    scalar q = Foam::pow(10.0, (12.0/nSections_)); 
    
    Info << "Volumetric spacing factor is " << q << nl << endl;
    
    Info << nl << "section" << token::TAB << "diameter" << token::TAB << "surface area" << token::TAB << "volume" << endl;
    
    for (label section = 0; section < d_.size(); section++)
    {
        v_[section] = v_[0]*Foam::pow(q,(section)); // calculating volume of all larger nodes using a geometric factor of 2. Note that volume of node 0 is the volume of a molecule.
        d_[section] = Foam::pow((6.0*v_[section]/M_PI),1.0/3.0);
        a_[section] = M_PI*Foam::pow(d_[section], 2.0);
        
        Info << section << token::TAB << d_[section] << token::TAB << a_[section] << token::TAB << v_[section] << endl;
    }
    
    Info << endl;
}


template<class CloudType>
const Foam::wordList Foam::SectionalDropletSynthesis<CloudType>::psdPropertyNames() const
{
    wordList sectionNames;
    sectionNames.setSize(nSections_);
    
    forAll(sectionNames,section)
    {
        word sectionName;
        
        if (section < 10)
        {
            sectionNames[section] = "section" + Foam::name(0) + Foam::name(section);
        }
        else
        {
            sectionNames[section] = "section" + Foam::name(section);
        }
    }
        
    return sectionNames;
}


template<class CloudType>
const Foam::wordList Foam::SectionalDropletSynthesis<CloudType>::physicalAerosolPropertyNames() const
{
    wordList physicalAerosolPropertyNames;
    physicalAerosolPropertyNames.setSize(nPhysicalAerosolProperties_);
    
    physicalAerosolPropertyNames[0] = "saturation";
    physicalAerosolPropertyNames[1] = "nucleationRate";
    physicalAerosolPropertyNames[2] = "condensationRate";
        
    return physicalAerosolPropertyNames;
}


template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::setSizeSplittingOperator()
{
    label nSecSqr = nSections_*nSections_;
    
    X_.setSize(Foam::pow(nSections_,label(3)), 0.0);
    
    for (label k=1; k <= nSections_-1; k++)
    {
        for (label i=0; i <= nSections_-1; i++)
        {
            for (label j=0; j <= nSections_-1; j++)
            {
                //Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
                if (v_[k] <= (v_[i]+v_[j]) && (v_[i]+v_[j]) < v_[k+1])
                {
                    X_[i*nSecSqr+j*nSections_+k] = (v_[k+1]-v_[i]-v_[j])/(v_[k+1]-v_[k]);
                }             
                else
                {
                    if (v_[k-1] <= (v_[i]+v_[j]) && (v_[i]+v_[j]) < v_[k])
                    {
                        X_[i*nSecSqr+j*nSections_+k] = (v_[i]+v_[j]-v_[k-1])/(v_[k]-v_[k-1]);                      
                    }
                    else 
                    {
                        X_[i*nSecSqr+j*nSections_+k] = 0;
                    }
                }              
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::particleInception(particleType& particle) 
{
    CloudType& cloud(this->owner());
    scalar deltaT = cloud.mesh().time().deltaT().value();
    
    scalar Ps = saturationPressure(particle.T());
    scalar ns = Ps/(kb_*particle.T());
    
    scalar S = 1.0;
        
    if (cloud.mesh().time().timeIndex() == 1)
    {
        S = 1.001;
        
        particle.psdProperties()[0] = S*ns;
    }
    
    S = particle.psdProperties()[0]/ns; // saturation ratio.

    scalar sigma = surfaceTension(particle.T()); 
    
    scalar vstar = criticalVolume(sigma, particle.T(), S);
    
    scalar kstar = criticalNode(sigma, particle.T(), S);
    scalar Jk = nucleationRate(sigma, particle.T(), S, ns);
    
    scalarList zeta(nSections_, 0.0);
        
    // operator to put nucleated particles in the bin just higher than k*.
    for (label k=1; k <= nSections_-1; k++)
    {
        if (vstar < v_[0])
        {
            zeta[1] = vstar/v_[1]; 
        }
        else if ((v_[k-1] <= vstar) && (vstar < v_[k])) 
        {
            zeta[k] = vstar/v_[k]; 
        }
        else
        {
            zeta[k] = 0.0;
        }
    }
    
    for(label k=1; k <= nSections_-2; k++)
    {
        particle.psdProperties()[k] += deltaT*Jk*zeta[k]; // changing PSD due to nucleation
    }
    
    particle.psdProperties()[0] -= deltaT*Jk*kstar;
}


template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::particleSurfaceGrowth(particleType& particle) 
{
    CloudType& cloud(this->owner());
    scalar deltaT = cloud.mesh().time().deltaT().value();
    
    scalar beta[nSections_][nSections_];
    
    scalar Ps = saturationPressure(particle.T());
    scalar ns = Ps/(kb_*particle.T());
    
    scalar S = 1.0;
        
    if (cloud.mesh().time().timeIndex() == 1)
    {
        S = 1.001;
        
        particle.psdProperties()[0] = S*ns;
    }
    
    S = particle.psdProperties()[0]/ns; // saturation ratio.

    scalar sigma = surfaceTension(particle.T()); 
    
    for(label j=0; j <= nSections_-1; j++)
    {
        scalar temp1 = 1.0/v_[0] + 1.0/v_[j];
        scalar temp2 = Foam::pow(v_[0],1.0/3.0) + Foam::pow(v_[j],1.0/3.0);
         
        beta[0][j] = Foam::pow(3.0/(4.0*M_PI), 1.0/6.0)*Foam::pow((6.0*kb_*particle.T()/liquidDensity_),0.5)*Foam::pow(temp1,0.5)*Foam::pow(temp2,2.0);
    }
    
    scalar molecularWeight = speciesThermo_[condensingSpeciesIndex_].W();
    
    scalarList N1s(nSections_, 0.0);
    
    for(label i=0; i <= nSections_-1; i++)
    {
        N1s[i] = ns*Foam::exp(4.0*sigma*molecularWeight/(Ru_*liquidDensity_*particle.T()*d_[i]));
    }
    
    scalar t1 = 0;
    scalar t2 = 0;
    scalar t3 = 0;
    scalar t4 = 0;
    
    for(label k=1; k <= nSections_-2; k++)
    {
        scalar addterm = 0.0;
        scalar subterm = 0.0;
        
        if (particle.psdProperties()[0] > N1s[k-1])
        {
            if (k==1)
            {
                addterm = 0.0;
            }
            else
            {
                addterm = (v_[0]/(v_[k]-v_[k-1]))*beta[0][k-1]*(particle.psdProperties()[0]-N1s[k-1])*particle.psdProperties()[k-1]; // growth of k due to condensation of monomers on k-1
                t1 += beta[0][k-1]*(particle.psdProperties()[0]-N1s[k-1])*particle.psdProperties()[k-1]; // loss of monomers that have condensed.
            }
        }        
        
        if (particle.psdProperties()[0] < N1s[k+1])
        {
            addterm = -(v_[0]/(v_[k+1]-v_[k]))*beta[0][k+1]*(particle.psdProperties()[0]-N1s[k+1])*particle.psdProperties()[k+1]; // growth of k due to evaporation of monomers from k+1
            t2 += beta[0][k+1]*(-particle.psdProperties()[0] + N1s[k+1])*particle.psdProperties()[k+1]; // gain of monomers that have evaporated.
        }
        
        if (particle.psdProperties()[0] < N1s[k])
        {
            subterm = -(v_[0]/(v_[k]-v_[k-1]))*beta[0][k]*(particle.psdProperties()[0]-N1s[k])*particle.psdProperties()[k]; // loss of k due to evaporation of monomers from k
            t3 += beta[0][k]*(-particle.psdProperties()[0]+N1s[k])*particle.psdProperties()[k]; // gain of monomers that have evaporated.
        }
        
        if (particle.psdProperties()[0] > N1s[k])
        {
            subterm = (v_[0]/(v_[k+1]-v_[k]))*beta[0][k]*(particle.psdProperties()[0]-N1s[k])*particle.psdProperties()[k]; // loss of k due to condensation of monomers on k          
            t4 += beta[0][k]*(particle.psdProperties()[0]-N1s[k])*particle.psdProperties()[k]; // loss of monomers that have condensed.
        }
        
        particle.psdProperties()[k] += deltaT*(addterm - subterm); // changing particle size distribution due to surface growth/evaporation
    }
    
    particle.psdProperties()[0] -= deltaT*(t1 + t4 - t3 - t2);
}


template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::particleCoagulation(particleType& particle) 
{
    CloudType& cloud(this->owner());
    scalar deltaT = cloud.mesh().time().deltaT().value();
    
    label nSecSqr = nSections_*nSections_;
    
    scalar beta[nSections_][nSections_];
    
    if (coagulationModel_ == "freeMolecularCollisionKernel")
    {
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                scalar v1 = 1.0/v_[i] + 1.0/v_[j];
                scalar v2 = Foam::pow(v_[i], 1.0/3.0) + pow(v_[j], 1.0/3.0);
                
                beta[i][j] = Foam::pow(3.0/(4.0*M_PI), 1.0/6.0)*Foam::pow((6.0*kb_*particle.T()/liquidDensity_), 0.5)*pow(v1, 0.5)*pow(v2, 2.0);               
            }
        }
    }
    else if (coagulationModel_ == "FuchsCollisionKernel")
    {
        scalar mu = dynamicViscosity(particle.T());
        
        scalar lambda = meanFreePath(particle.T(), particle.pc(), mu);
        
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                scalar Kn1 = (2.0*lambda)/d_[i];
                scalar Kn2 = (2.0*lambda)/d_[j];
                
                scalar D1 = (kb_*particle.T())/(3.0*M_PI*mu*d_[i])*((5.0 + 4.0*Kn1 + 6.0*Kn1*Kn1 + 18.0*Kn1*Kn1*Kn1)/(5.0 - Kn1 + (8.0 + M_PI)*Kn1*Kn1));
                scalar D2 = (kb_*particle.T())/(3.0*M_PI*mu*d_[j])*((5.0 + 4.0*Kn2 + 6.0*Kn2*Kn2 + 18.0*Kn2*Kn2*Kn2)/(5.0 - Kn2 + (8.0 + M_PI)*Kn2*Kn2));
                
                scalar mp1 = v_[i]*liquidDensity_;
                scalar mp2 = v_[j]*liquidDensity_;
                
                scalar c1 = sqrt((8.0*kb_*particle.T())/(M_PI*mp1));
                scalar c2 = sqrt((8.0*kb_*particle.T())/(M_PI*mp2));
                
                scalar l1 = (8.0*D1)/(M_PI*c1);
                scalar l2 = (8.0*D2)/(M_PI*c2);
                
                scalar g1 = (Foam::pow((d_[i] + l1), 3) - pow((d_[i]*d_[i] + l1*l1), 1.5))/(3.0*d_[i]*l1) - d_[i];
                scalar g2 = (Foam::pow((d_[j] + l2), 3) - pow((d_[j]*d_[j] + l2*l2), 1.5))/(3.0*d_[j]*l2) - d_[j];
                
                beta[i][j] = 2.0*M_PI*(D1 + D2)*(d_[i] + d_[j])/((d_[i] + d_[j])/(d_[i] + d_[j] + 2.0*sqrt(g1*g1 + g2*g2)) + (8.0*(D1 + D2))/(sqrt(c1*c1 + c2*c2)*(d_[i] + d_[j])));
            }
        }
    }
    else
    {
        for (label i=0; i <= nSections_-1; i++)
        {
            for(label j=0; j <= nSections_-1; j++)
            {
                beta[i][j] = 0.0;               
            }
        }
    }  
    
    for(label k=1; k <= nSections_-2; k++)
    {
        scalar add = 0.0; // Addition term when i and j collide to form a k sized particle.
        scalar sub = 0.0; // Subtraction term when k collides with any other particle.
        
        for (label i=1; i <= nSections_-1; i++)
        {
            sub += beta[k][i]*particle.psdProperties()[i];
            
            for (label j=1; j<=k; j++)
            {
                add += X_[i*nSecSqr+j*nSections_+k]*beta[i][j]*particle.psdProperties()[i]*particle.psdProperties()[j];
            }
        }
        
        particle.psdProperties()[k] += deltaT*(0.5*add - particle.psdProperties()[k]*sub);
    }
}


template <class CloudType>
void Foam::SectionalDropletSynthesis<CloudType>::gasPhaseCooling(particleType& particle) 
{
    CloudType& cloud(this->owner());
    scalar deltaT = cloud.mesh().time().deltaT().value();
    
    particle.hA() = particle.hA() - 520.4*coolrate_*deltaT;
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::surfaceTension(scalar T) 
{
    scalar A_ = 948.0;
    scalar B_ = 0.202;
    
    return (A_ - B_*T)*1e-3;
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::saturationPressure(scalar T) 
{
    scalar C_ = 13.07;
    scalar D_ = 36373.0;
    
    return Foam::exp(C_ - D_/T)*101325.0;
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::nucleationRate(scalar sigma, scalar T, scalar S, scalar ns) 
{
    scalar s1 = M_PI*Foam::pow(6*v_[0]/M_PI,2.0/3.0);
    scalar m1 = liquidDensity_*v_[0];
    
    scalar theta = (s1*sigma)/(kb_*T);
    
    scalar a = (2*sigma)/(M_PI*m1); 
    scalar b = theta - (4*Foam::pow(theta,3.0))/(27.0*Foam::pow(Foam::log(S),2.0)); 
    
    return Foam::pow(ns, 2.0)*S*v_[0]*Foam::pow(a, 0.5)*Foam::exp(b); 
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::criticalNode(scalar sigma, scalar T, scalar S) 
{
    scalar s1 = M_PI*Foam::pow(6*v_[0]/M_PI,2.0/3.0);
        
    scalar theta = (s1*sigma)/(kb_*T);
    
    scalar c = 2.0/3.0*theta/Foam::log(S); 
    
    return Foam::pow(c, 3.0); 
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::criticalVolume(scalar sigma, scalar T, scalar S) 
{
    scalar dpstar = 4.0*sigma*v_[0]/(kb_*T*Foam::log(S)); 
    
    return M_PI*Foam::pow(dpstar, 3.0)/6.0;
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::dynamicViscosity(scalar T) 
{
    // calculating dynamic viscosity following the Sutherland law
    scalar C1 = 1.966e-06;
    scalar S  = 147.47;
        
    return C1*Foam::pow(T, 1.5)/(T + S);
}


template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::meanFreePath(scalar T, scalar p, scalar mu) 
{
    return mu/p*Foam::sqrt(M_PI*Ru_/1000*T/(2.0*0.04));
}
        

template <class CloudType>
Foam::scalar Foam::SectionalDropletSynthesis<CloudType>::calcMolecularWeightMixture(scalarField& Y) 
{
    scalar molecularWeightMixture = 0.0;
    
    forAll(speciesThermo_,specieI) 
    {
        molecularWeightMixture += Y[specieI]/speciesThermo_[specieI].W();
    }
        
    return 1.0/molecularWeightMixture;    
}

// ************************************************************************* //
