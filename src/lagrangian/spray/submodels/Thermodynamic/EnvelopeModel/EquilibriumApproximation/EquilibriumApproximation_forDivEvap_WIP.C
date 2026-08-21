/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "EquilibriumApproximation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EquilibriumApproximation<CloudType>::EquilibriumApproximation
(
    const dictionary& dict,
    CloudType& owner
)
:
    EnvelopeModel<CloudType>(dict,owner,typeName)
{
    this->setSwitches();
    setTableInfo();
}


template<class CloudType>
Foam::EquilibriumApproximation<CloudType>::EquilibriumApproximation
(
    const EquilibriumApproximation<CloudType>& ea
)
:
    EnvelopeModel<CloudType>(ea)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EquilibriumApproximation<CloudType>::~EquilibriumApproximation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class CloudType>
inline Foam::scalar& Foam::EnvelopeModel<CloudType>::nOxNet()
{
    return nOxNet_;
}

template<class CloudType>
inline Foam::labelField& Foam::EnvelopeModel<CloudType>::n0Fuel()
{
    return n0Fuel_;
}

template<class CloudType>
inline Foam::labelField& Foam::EnvelopeModel<CloudType>::n1()
{
    return n1_;
}
*/
template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::setTableInfo()
{
    const dictionary elementsDict(this->coeffDict().subDict("elements"));
    const dictionary inertFuelSpeciesDict(this->coeffDict().subDict("inertFuelSpecies"));
    const dictionary speciesElementsDict(this->coeffDict().subDict("speciesElements"));
    const word oxidizer(this->coeffDict().lookup("oxidizer"));
    const dictionary speciesStDict(this->coeffDict().subDict("stoichiometry"));

    printOut_ = this->coeffDict().lookupOrDefault("printOut", false);
    oxidizer_ = oxidizer; 
    oxRatio_ = 0.0;
    nElements_ = elementsDict.size();
    nSpeciesSt_ = speciesStDict.size() + inertFuelSpeciesDict.size();
    inertStart_ = speciesStDict.size();
    fuelRatio_.setSize(nElements_,0.0);
    stThermoId_.setSize(nSpeciesSt_,-1);
    speciesSt_.setSize(nSpeciesSt_);

    elementNames_ = elementsDict.toc();
    inertSpecies_.setSize(this->owner().pope().composition().componentNames().size(),false);
    elements_.setSize(nElements_);
    forAll(elements_,ne)
    {
        elements_[ne].setSize(this->owner().pope().composition().componentNames().size());
        if(elementsDict.found(elementNames_[ne]))
        {
            fuelRatio_[ne] = readScalar(elementsDict.lookup(elementNames_[ne]));
        }
    }

    label countSt = 0;
    forAll(this->owner().pope().composition().componentNames(), ns)
    {
        const word currentSpecie = this->owner().pope().composition().componentNames()[ns];
        if(inertFuelSpeciesDict.found(currentSpecie))
        {
            inertSpecies_[ns] = true;

        }

        if(speciesElementsDict.found(currentSpecie))
        {
            List<scalar> currentElements(speciesElementsDict.lookup(currentSpecie));
            forAll(currentElements, ne)
            {
                elements_[ne][ns] =  currentElements[ne];
                if(currentSpecie == oxidizer_) oxRatio_ += currentElements[ne];
            }
        }
	    else
        {
            forAll(elements_,ne)
            {
                elements_[ne][ns] = 0;
            }
        }

	    if(speciesStDict.found(currentSpecie) || inertFuelSpeciesDict.found(currentSpecie))
        {
            stThermoId_[countSt] = ns;
	        speciesSt_[countSt] = currentSpecie;
	        countSt++;
        }
/*
	if(inertFuelSpeciesDict.found(currentSpecie))
        {
            stThermoId_[countSt] = ns;
            speciesSt_[countSt] = currentSpecie;
            countSt++;
        }
*/
    }

    Info << "printOut? " << printOut_ << endl;
    Info << "nElements_: " << nElements_ << ", elements_.size(): " << elements_.size() << endl;
    Info << "inertSpecies: " << inertSpecies_ << endl;
    Info << nSpeciesSt_ << ": stoichiometrix species: " << speciesSt_ << endl;
    Info << "stoichiometrix IDs " << stThermoId_ << endl;
    Info << "elements_:\n" << elements_ << endl;
    Info << "oxidizer: " << oxidizer_ << endl;
    Info << "fuelRatio: " << fuelRatio_ << ", oxRatio: " << oxRatio_ << endl;
}

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::initialize
(
    const scalarField Y0,
    const scalar T0,
    const label sCell
)
{
    envelopeSwitch_ = this->checkState(T0);
    if(envelopeSwitch_) this->gasTable(Y0);
    if(nOxNet_ <= 0.0) envelopeSwitch_ = false;
    fSurfOld_ = 1.0;
    Info << "envelopeSwitch: " << envelopeSwitch_ << endl;
}

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::gasTable
(
    const scalarField Y0
)
{
    // Initialise variables
    scalar nSpecies = this->owner().pope().composition().componentNames().size();
    MW0_ = this->owner().pope().composition().mixtureMW(Y0);
    scalarField X0 = this->owner().pope().composition().X(Y0);

    // In 1 mol of gas (sp), the quantity of each element and O2
    scalar nOx = 0.0;
    scalarField n0Fuel;
    n0_.setSize(nElements_);
    n0Fuel.setSize(nElements_,0.0);

    n0_ = Zero;
    // Populating each group
    forAll(X0,ns)
    {
        if(this->owner().pope().composition().componentNames()[ns] == oxidizer_) nOx = X0[ns];

        forAll(elements_,ne)
        {
            n0_[ne] += elements_[ne][ns] * X0[ns];
            if(!inertSpecies_[ns]) n0Fuel[ne] += elements_[ne][ns] * X0[ns];
//	    else n0Fuel[ne] += elements_[ne][ns] * X0[ns];
	    }
    }

    // In 1 mol of gas (sp), additional O-element required to convert all radicals into its product
    scalar nOxLoss = 0.0;
    forAll(elements_,ne)
    {
        nOxLoss += n0Fuel[ne] * fuelRatio_[ne];
    }
    nOxLoss /= oxRatio_;
    // In 1 mol of gas (sp), O2 left after the conversion 
    nOxNet_ = nOx - nOxLoss;
}

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::solveTable
(
    const scalarField Y1
)
{
    scalar MW1 = this->owner().pope().composition().mixtureMW(Y1);
    scalarField X1 = this->owner().pope().composition().X(Y1);
    scalarField n1;
    n1.setSize(nElements_,0.0);
    YSt_.setSize(this->owner().pope().composition().componentNames().size(),0.0);
    label nElementsOx = 0;
    forAll(X1,ns)
    {
        forAll(elements_,ne)
        {
            n1[ne] += elements_[ne][ns] * X1[ns];
//	    Info << this->owner().pope().composition().componentNames()[ns] << " " << oxidizer_ << endl;
            if(this->owner().pope().composition().componentNames()[ns] == oxidizer_) nElementsOx += elements_[ne][ns];
        }

//        if(!inertSpecies_[ns]) fSt_ += Y1[ns];
    }

    scalar nOxReq = 0.0;
    forAll(elements_,ne)
    {
        nOxReq += n1[ne] * fuelRatio_[ne];
    }
    if(printOut_) Info << "nOxReq: " << nOxReq << ", nOxNet: " << nOxNet_ << endl;
    scalar n0Req = nOxReq / nOxNet_;
    if(printOut_) Info << "n0Req=nOxReq/nOxNet: " << n0Req << "=" << nOxReq << "/" << nOxNet_ << endl;

    scalarField nProd;
    label nOxElementsProd = 0;
//    label nElementsProd = 0;
    nProd.setSize(nSpeciesSt_,0.0);
    scalar mProd = 0.0;

    forAll(nProd,np)
    {
        if(printOut_) Info << "( ";
        forAll(elements_,ne)
        {
            if(fuelRatio_[ne]>=0.0 && np < inertStart_)
            {
                nProd[np] += elements_[ne][stThermoId_[np]] * n1[ne];
//		nElementsProd += elements_[ne][stThermoId_[np]];
                if(printOut_) Info << " +(" << elements_[ne][stThermoId_[np]] <<"*" << n1[ne] << ") [" << nProd[np] << "] " << elementNames_[ne];
            }
	    else if(np > inertStart_-1)
            {
                if(fuelRatio_[ne]==0.0) 
                {
                    nProd[np] += n0Req * n0_[ne];
                    if(printOut_) Info << " +(" << n0Req << "*" << n0_[ne]  << ") [" << nProd[np] << "] "  << elementNames_[ne];
                }
		else
		{
                    if(printOut_) Info << " +(0*0) [" << nProd[np] << "] "  << elementNames_[ne];
		}
            }
	    else
            {
                nProd[np] += n0Req * n0_[ne];
                nOxElementsProd = elements_[ne][stThermoId_[np]];
                if(printOut_) Info << " +(" << n0Req << "*" << n0_[ne] << ") [" << nProd[np] << "] "  << elementNames_[ne];
            }
//            if(printOut_) Info << "nOxElementsProd: " << nOxElementsProd << ", nElementsOx: " << nElementsOx << endl;
        }
//        if(printOut_) Info << " nProd: " << nProd[np] << endl;
	if(nElementsOx>0.0 && np < inertStart_) 
	{
	    nProd[np] = nProd[np] * nOxElementsProd / nElementsOx;
            if(printOut_) Info << ") * " << nOxElementsProd << "/" << nElementsOx << " nProd: " << nProd[np] <<endl;
	}
	if(nElementsOx>0.0 && np > inertStart_-1)
	{
	    nProd[np] = nProd[np] * 1 / nElementsOx;
            if(printOut_) Info << ") * " << "1/" << nElementsOx << " nProd: " << nProd[np] <<endl;
	}
//        if(nElementsOx>0.0) nProd[np] = nProd[np] / nElementsProd;
//exclude oxidizer elements        nProd[np] = elements_[ne][stThermoId_[np]];
        YSt_[stThermoId_[np]] = nProd[np] * this->owner().pope().composition().molWt(stThermoId_[np]);
	mProd += YSt_[stThermoId_[np]];
    }
    if(printOut_) Info << "mProd: " << mProd << endl;
    fSt_ = 1.0 * MW1 / (1.0 * MW1 + n0Req * MW0_);
    YSt_ /= mProd;
    if(printOut_) Info << "fSt:  " << fSt_ << endl; //"\n YSt: " << YSt_ << endl;


//////////////////////////////////////////////////////////////////////////////////////
// To consume 1 mol of volatile, quantity of pure O2 required
/*   
    scalar n_o2_req = 0.5 * (n_c_vol * 2.0 + n_h_vol * 0.5 - n_o_vol);

    // To consume 1 mol of volatile, quantity of gas (sp) required
    scalar n_sp_req = n_o2_req / n_o2_net;
    // and its product
// <<<<< currenttly here
    scalar n_co2_prod = n_c_vol + n_sp_req * n_c_sp;
    scalar n_h2o_prod = 0.5 * (n_h_vol + n_sp_req * n_h_sp);
    scalar n_n2_prod = 0.5 * (n_n_vol + n_sp_req * n_n_sp);

    // Balanced reaction: 
    // 1 {VOL} + n_sp_req {SP} ==> n_co2_prod {CO2} + n_h2o_prod {H2O} + n_n2_prod {N2}
    // Compute local stoich. mixture fraction and local stoich. composition (final product)
    fst_ = 1.0 * MW_vol / (1.0 * MW_vol + n:q
    _sp_req * MW_sp);
    scalar mp = n_co2_prod * 44.0 + n_h2o_prod * 18.0 + n_n2_prod * 28.0;
    scalar y_co2_st = n_co2_prod * 44.0 / mp;
    scalar y_h2o_st = n_h2o_prod * 18.0 / mp;
    scalar y_n2_st = 1.0 - y_co2_st - y_h2o_st;

    // Populate the equilibrium envelope table at stoichiometry
    scalarField Y_stoich(nSpecies, 0.0);
    forAll(td.cloud().pope().composition().componentNames(), ns)
    {
        if (td.cloud().pope().composition().componentNames()[ns] == "CO2")
        {
            Y_stoich[ns] = y_co2_st;
        }
        else if (td.cloud().pope().composition().componentNames()[ns] == "H2O")
        {
            Y_stoich[ns] = y_h2o_st;
        }
        else if (td.cloud().pope().composition().componentNames()[ns] == "N2")
        {
            Y_stoich[ns] = y_n2_st;
        }
    }
    
    // Transfer the table into its final form.
    Yst_ = Y_stoich;    // f = fst_
    f0_ = fGas_;
    f1_ = 1.0;
*/
}

template<class CloudType>
//**//Foam::scalarField Foam::EquilibriumApproximation<CloudType>::fSurface
Foam::scalar Foam::EquilibriumApproximation<CloudType>::fSurface //**//
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar f0,
    const scalarField YSurf
)
{
//**//    scalarField fSurf;
//**//    fSurf.setSize(this->owner().linkFG().size(),-2.0);
    scalar fSurf = -2.0; //**//
//  if(Y1Liq != Y1LiqOld) this->solveTable(this->Y1(Y1Liq));
    this->solveTable(this->Y1(Y1Liq));
    if(envelopeSwitch_ && fSurfOld_ >= fSt_) fSurf = this->fNoEnvelope(YSt_,Y1Liq,fSt_,YSurf);
    else fSurf = this->fNoEnvelope(Y0,Y1Liq,f0,YSurf);
    if(printOut_) Info << "fSurf: " << fSurf << endl;
//**//    fSurfOld_ = fSurf[this->owner().linkFG().size()-1];
    fSurfOld_ = fSurf;
    return fSurf;
}




template<class CloudType>
Foam::scalarField Foam::EquilibriumApproximation<CloudType>::YSurface
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar f0,
    const scalar fSurf
) const
{
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YSurf;
    YSurf.setSize(Y1.size(),0.0);
    if(envelopeSwitch_ && fSurf >= fSt_) YSurf = this->YNoEnvelope(YSt_,Y1,fSt_,fSurf);
    else YSurf = this->YNoEnvelope(Y0,Y1,f0,fSurf);
//    if(printOut_) Info << "YSurf: " << YSurf << endl;
    return YSurf;
}

template<class CloudType>
Foam::scalarField Foam::EquilibriumApproximation<CloudType>::YFilm
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar f0,
    const scalar fFilm
) const
{   
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YFilm;
    YFilm.setSize(Y1.size(),0.0);
    if(printOut_) Info << "fFilm: " << fFilm << endl;
    if(printOut_) Info << "YSt: " << YSt_ << endl;
    if(printOut_) Info << "Y0: " << Y0 << endl;
    if(envelopeSwitch_ && fSurfOld_ > fSt_ && fFilm >= fSt_) YFilm = this->YNoEnvelope(YSt_,Y1,fSt_,fFilm);
    else if(envelopeSwitch_ && fSurfOld_ > fSt_ && fFilm < fSt_) YFilm = this->YEnvelope(Y0,YSt_,f0,fSt_,fFilm);
    else YFilm = this->YNoEnvelope(Y0,Y1,f0,fFilm);
    if(printOut_) Info << "YFilm: " << YFilm << endl;
    return YFilm;
}

template<class CloudType>
Foam::scalarField Foam::EquilibriumApproximation<CloudType>::YEnvelope
(
     const scalarField Y0,
     const scalarField Y1,
     const scalar f0,
     const scalar f1,
     const scalar fX
) const
{
    scalarField YX;
    YX.setSize(Y0.size(),0.0);
    YX = Zero;

    forAll(YX, ns)
    {
        YX[ns] = (Y1[ns] - Y0[ns]) / (f1 - f0 + VSMALL) * (fX - f0) + Y0[ns];
//	if(printOut_) Info << "(" << Y1[ns] << " - " <<  Y0[ns] << ") / (" << f1 << " - " << f0 << " + " << VSMALL << ") * (" << fX << " - " << f0 << ") + " << Y0[ns] << endl;
    }
    return YX;
}

template<class CloudType>
Foam::scalar Foam::EquilibriumApproximation<CloudType>::haSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar ha0,
    const scalar ha1,
    const scalar fSurf
) const
{
    scalar haSurf = -2.0;
    if(envelopeSwitch_ && fSurf >= fSt_) haSurf = this->haEnvelope(ha0,ha1,fSurf);
    else haSurf = this->haNoEnvelope(YSurf,p,TD);

    return haSurf;
}

template<class CloudType>
Foam::scalar Foam::EquilibriumApproximation<CloudType>::haEnvelope
(
     const scalar ha0,
     const scalar ha1,
     const scalar fSurf
) const
{
    scalar haX;

//**//    scalar ratio = (fSurf - fSt_) / (this->f1_[this->owner().linkFG().size()-1] - fSt_ + ROOTVSMALL);
    scalar ratio = (fSurf - fSt_) / (this->f1_ - fSt_ + ROOTVSMALL); //**//
    haX = ha0 - ratio * (ha0 - ha1);
    return haX;
}

template<class CloudType>
Foam::scalar Foam::EquilibriumApproximation<CloudType>::TEnvelope
(
     const scalarField YSurf,
     const scalar haSurf,
     const scalar p,
     const scalar T0
) const
{
    scalar TSurf;

    TSurf = this->owner().pope().composition().particleMixture(YSurf).THa(haSurf, p, T0);

    return TSurf;
}

//-Surface temperature
template<class CloudType>
Foam::scalar Foam::EquilibriumApproximation<CloudType>::TSurface
(
    const scalarField YSurf,
    const scalar p,
    const scalar TD,
    const scalar TGas,
    const scalar fSurf,
    const scalar haSurf
) const
{
    scalar TSurf = -2.0;
    if(envelopeSwitch_ && fSurf >= fSt_) TSurf = this->TEnvelope(YSurf, haSurf, p, TD);
    else TSurf = this->TNoEnvelope(TD);

    return TSurf;
}

// ************************************************************************* //
