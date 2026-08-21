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

#include "PremixedMixingPopeCloud.H"
#include "CloudMixingModel.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::setModels(const mmcVarSet& Xi)
{
    mixingModel_.reset
    (
        CloudMixingModel<PremixedMixingPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this,
            Xi //Since reference variables are inside submodel!!!
        ).ptr()
    );
}


template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::cloudReset(PremixedMixingPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
    
    mixingModel_.reset(c.mixingModel_.ptr());
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PremixedMixingPopeCloud<CloudType>::PremixedMixingPopeCloud
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
  
    mixingPopeCloud(),
    
    cloudCopyPtr_(nullptr),
    
    mixingModel_(nullptr),
    
    XiRrelaxationDict_(this->cloudProperties().subDict("XiRrelaxation")),
    
    XiRrelaxEnabled_(XiRrelaxationDict_.lookup("enabled")),

    XiRrelaxType_(XiRrelaxationDict_.lookupOrDefault<word>("relaxationType", "dynamicRelax")), 

    XiRrelaxConstValue_(XiRrelaxationDict_.lookupOrDefault<scalar>("constRelFactor", 1.0)), 

    XiRDynRelaxLower_(XiRrelaxationDict_.lookupOrDefault<scalar>("cAverLowerLim", 0.475)), 

    XiRDynRelaxUpper_(XiRrelaxationDict_.lookupOrDefault<scalar>("cAverUpperLim", 0.525)), 

    XiRDynRelaxRate_(XiRrelaxationDict_.lookupOrDefault<scalar>("dynChangeRate", 0.1)), 

    XiRDynStepNum_(XiRrelaxationDict_.lookupOrDefault<label>("relaxStepNum", 1000)), 

    XiRDynRelaxMin_(XiRrelaxationDict_.lookupOrDefault<scalar>("dynRFmin", 0.01)), 

    XiRDyncPcrit_(XiRrelaxationDict_.lookupOrDefault<scalar>("cPcrit", 0.8)), 

    XiRDyncPcritMax_(XiRrelaxationDict_.lookupOrDefault<scalar>("cPcritMax", 0.9)),  

    XiRDyncPcritMin_(XiRrelaxationDict_.lookupOrDefault<scalar>("cPcritMin", 0.4))
    
{
    Info << "Creating premixed Mixing Pope Particle Cloud." << nl << endl;

    setModels(Xi); // passing mmcVarSet
    
    Info << nl << "Mixing model constructed." << endl;
    
    setEulerianStatistics();
    
    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Mixing Pope particle cloud data from file." << endl;
            
            particleType::mixingParticleIOType::readFields(*this, this->mixing());
        }
    
        else
        {
            if(initAtCnstr)
            {
                Info << "Initial realease of Pope particles into the finite volume field." << nl << endl;
            
                this->initReleaseParticles();
            }
        }
    }
    
    if (!XiRrelaxEnabled_)
    {
    	Info << nl << "Reference variable relaxation is disabled." << endl;
    }
    else
    {
    	Info << nl << "Reference variable relaxation is enabled." << endl;
        
        if (XiRrelaxType_ == "constValue")
        {
        	Info << "Relaxation type: constValue" << nl
             	 << token::TAB << "constRelFactor:     " << XiRrelaxConstValue_ << nl << endl;
        }
        else if (XiRrelaxType_ == "dynamicRelax")
        {
        	Info << "Relaxation type: dynamicRelax" << nl
             	 << token::TAB << "cAverLowerLim:     " << XiRDynRelaxLower_ << nl
             	 << token::TAB << "cAverUpperLim:     " << XiRDynRelaxUpper_ << nl
             	 << token::TAB << "dynChangeRate:     " << XiRDynRelaxRate_ << nl
             	 << token::TAB << "relaxStepNum:      " << XiRDynStepNum_ << nl
             	 << token::TAB << "dynRFmin:          " << XiRDynRelaxMin_ << nl
             	 << token::TAB << "cPcrit:            " << XiRDyncPcrit_ << nl
             	 << token::TAB << "cPcritMax:         " << XiRDyncPcritMax_ << nl
             	 << token::TAB << "cPcritMin:         " << XiRDyncPcritMin_ << nl << endl;
             		 
            if (Pstream::myProcNo() == 0) 
            {
            	dynamicRelaxOut_.reset
        		(   
            		new OFstream("dynamicRelaxOut_"+cloudName+".dat")
        		);
        		dynamicRelaxOut_() << "Text file for debugging." << nl;
            }
		}
        else
        {
        	FatalError << "Selected reference variable relaxation type does not exist."<<nl
                       << "Valid types are: constValue, dynamicRelax"
                       << exit(FatalError);
        }
            
	}

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PremixedMixingPopeCloud<CloudType>::~PremixedMixingPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::setParticleProperties
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
    
    particle.dx() = 0;
    
    label numXiR = this->mixing().numXiR();

    //- Extension for a set of variables    
    particle.dXiR().setSize(numXiR,0.0);
    
    if (!iniRls)
    {
        particle.XiR() = this->mixing().XiR0(patchI,patchFace);
    }
    else
    {
        particle.XiR() = this->mixing().XiR0(particle.cell());
    } 
    
    particle.mixTime() = 0;
    particle.mixExt() = 0;
    particle.progVarInt() = 0;
    particle.XiRrlxFactor() = 1.0;
    particle.rlxStep() = 0;

}


template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::setEulerianStatistics()
{
    
    if (this->eulerianStatsDict().found("dx"))
    {
        const dimensionSet dim = dimLength;
        
        this->eulerianStats().newProperty("dx",dim);
    }
    
    //- statistics of reference variables (mixing distances)
    forAll(this->mixing().XiRNames(),XiRI)
    {
        word refVarName = "d" + mixing().XiRNames()[XiRI];
        
        if (this->eulerianStatsDict().found(refVarName))
        {
            const dimensionSet dim = dimless;// they will be considered dimless 
            
            this->eulerianStats().newProperty(refVarName,dim);   
        }
    }
}


template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();

    forAllIters(*this, iter)
    {
        this->eulerianStats().findCell(iter().position());

        //- statistics of distance for reference variables 
        forAll(mixing().XiRNames(),XiRI)
        {
            const word& refVarName = mixing().XiRNames()[XiRI];
            const word& dRefVarName = "d" + refVarName;
            
            if (this->eulerianStatsDict().found(dRefVarName))
                this->eulerianStats().calculate(dRefVarName,iter().wt(),iter().dXiR(refVarName));
        }
            
        if (this->eulerianStatsDict().found("dx"))
            this->eulerianStats().calculate("dx",iter().wt(),iter().dx());
    }
}


template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::mixingParticleIOType::writeFields(*this, this->mixing());
    }
}

template<class CloudType>
template<class TrackCloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename particleType::trackingData& td
)
{
	if ((this->mixing().numXiR() > 1) && XiRrelaxEnabled_)
	{
		FatalError << "Reference variable relaxation is implemented only for a SINGLE reference variable." 
                   << exit(FatalError);
	}
	else 
	{ 
		this->calcXiRrelax();
	}
	
    CloudType::solve(cloud,td);
}

template<class CloudType>
void Foam::PremixedMixingPopeCloud<CloudType>::calcXiRrelax()
{
	if (!XiRrelaxEnabled_)
    {
        forAllIters(*this, iter)
        {
        	particleType& p = iter();
        	p.XiRrlxFactor() = 1.0;
        }
    }
    else
	{
		if (XiRrelaxType_ == "constValue")
        {
			forAllIters(*this, iter)
            {
            	particleType& p = iter();
            	p.XiRrlxFactor() = XiRrelaxConstValue_;
            }
        }
        else if (XiRrelaxType_ == "dynamicRelax") 
        { 
            scalar XiRrelaxOut = 1.0;
            label stepCountOut = 0;

            // Assume first particle already existed so it has a correct relaxation factor on it
            const particleType& firstParticle = this->begin()();
            XiRrelaxOut = firstParticle.XiRrlxFactor();
            stepCountOut = firstParticle.rlxStep();

			reduce(XiRrelaxOut, minOp<scalar>());
            reduce(stepCountOut, maxOp<label>());
			
			scalar pWtSum = 0.0;
            scalar cPxpWt = 0.0;
            scalar avercP = 0.0;

			forAllIters(*this, iter)
			{   
				particleType& p = iter();

				if ((p.progVarInt() >= XiRDynRelaxLower_) && (p.progVarInt() <= XiRDynRelaxUpper_))
				{
					cPxpWt += p.wt() * p.XiC()[0];
					pWtSum += p.wt();
				}
			} 

            reduce(pWtSum, sumOp<scalar>());
            reduce(cPxpWt, sumOp<scalar>());

            avercP = cPxpWt / (pWtSum + VSMALL);

			if ((avercP > XiRDyncPcrit_) && (stepCountOut >= XiRDynStepNum_))
            {
            	XiRrelaxOut = XiRrelaxOut * XiRDynRelaxRate_;
            	stepCountOut = 1;
            }
            else if ((avercP > XiRDyncPcritMax_))
            {
            	XiRrelaxOut = XiRDynRelaxMin_;
            	stepCountOut = 1;
            }
            else if ((avercP < XiRDyncPcritMin_) && (stepCountOut >= 10))
            {
            	XiRrelaxOut = 1.0;
            	stepCountOut = 1;
            }
            else  stepCountOut++;

			if (XiRrelaxOut < XiRDynRelaxMin_)
            {
            	XiRrelaxOut = XiRDynRelaxMin_;
            }

			forAllIters(*this, iter)
			{
            	particleType& p = iter();
            	p.XiRrlxFactor() = XiRrelaxOut;
            	p.rlxStep() = stepCountOut;
            }

            if(Pstream::myProcNo() == 0) 
            {
            	dynamicRelaxOut_()  << XiRrelaxOut << token::TAB << stepCountOut << token::TAB << avercP << endl;
            	Info << "\tRelaxation factor: " << XiRrelaxOut << token::TAB << "relaxation step: "<< stepCountOut << endl;
            }
        }
        else 
        {
            FatalError << "Selected reference variable relaxation type does not exist."<<nl
                       << "Valid types are: constValue, dynamicRelax"
                       << exit(FatalError);
        }
    }
}

// ************************************************************************* //
