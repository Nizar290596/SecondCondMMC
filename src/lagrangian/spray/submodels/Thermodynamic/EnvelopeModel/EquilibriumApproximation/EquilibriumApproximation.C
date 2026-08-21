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

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::setTableInfo()
{
    const dictionary elementsDict(this->coeffDict().subDict("elements"));
    const dictionary inertFuelSpeciesDict(this->coeffDict().subDict("inertFuelSpecies"));
    const dictionary speciesElementsDict(this->coeffDict().subDict("speciesElements"));
    const word oxidizer(this->coeffDict().lookup("oxidizer"));
    const dictionary speciesStDict(this->coeffDict().subDict("stoichiometry"));

    printOut_ = this->coeffDict().lookupOrDefault("printOut", false);
    oxidizer_ = oxidizer;                //- name of oxidizer, usually O2
    oxRatio_ = 0.0;                      //- O2 consists of 2 O
    nElements_ = elementsDict.size();    //- number of considered elements, eg. 4: C, H, O, N
    nSpeciesSt_ = speciesStDict.size() + inertFuelSpeciesDict.size(); //- number of species at stoichiometry, eg. 3: CO2, H2O and N2
//    inertStart_ = speciesStDict.size();
    fuelRatio_.setSize(nElements_,0.0);  //- number of elemnts in product, in the current implementation every element (except O) can only appear in one product
    stThermoId_.setSize(nSpeciesSt_,-1); //- species IDs of products
    speciesSt_.setSize(nSpeciesSt_);     //- species names of products

    elementNames_ = elementsDict.toc();  //- names of the elements
//    inertSpecies_.setSize(this->owner().pope().composition().componentNames().size(),false);
    stSpecies_.setSize(this->owner().pope().composition().componentNames().size(),false); //- boolean list of species at stoichiometry
    elements_.setSize(nElements_);       //- array with number of considered elements for each species, maybe this data can be accessed from the chemkin reader

    //- read and construct all the necessary data
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
//        if(inertFuelSpeciesDict.found(currentSpecie))
//        {
//            inertSpecies_[ns] = true;
//        }

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
            stSpecies_[ns] = true;
	        countSt++;
        }
    }

    //- print information about the read data
    Info << "printOut? " << printOut_ << endl;
    Info << "nElements_: " << nElements_ << ", elements_.size(): " << elements_.size() << endl;
//    Info << "inertSpecies: " << inertSpecies_ << endl;
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
    envelopeSwitch_ = this->checkState(T0);     //- check if envelope conditions are satisfied
    if(envelopeSwitch_) this->gasTable(Y0);     //- generate gas side of envelope table
    if(nOxNet_ <= 0.0) envelopeSwitch_ = false; //- check if there's enough oxidizer for an envelope flame
    fSurfOld_ = 2.0;
    fSt_ = -1.0;
    Y1LiqOld_.setSize(this->owner().linkFG().size(),-1.0);
    forAll(Y1LiqOld_,nfs) Y1LiqOld_[nfs]=-1.0;  //- IS THIS NEEDED? It SHOULDN'T
    if(printOut_) Info << "envelopeSwitch: " << envelopeSwitch_ << endl;
}

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::gasTable
(
    const scalarField Y0
)
{
    // Initialise variables
    MW0_ = this->owner().pope().composition().mixtureMW(Y0);
    scalarField X0 = this->owner().pope().composition().X(Y0);

    // In 1 mol of gas (sp), the quantity of each element and oxidizer
    scalar nOx = 0.0;               //- mole fraction of oxidizer
    scalarField n0Fuel;             //- number of fuel elements in sp
    n0_.setSize(nElements_,0.0);    //- number of elements in sp
    n0Fuel.setSize(nElements_,0.0);

    //- Populating each group (nOx and n0Fuel)
    forAll(X0,ns)
    {
        const word currentSpecie = this->owner().pope().composition().componentNames()[ns];
        if(currentSpecie == oxidizer_) nOx = X0[ns];

        forAll(elements_,ne)
        {
            n0_[ne] += elements_[ne][ns] * X0[ns];
//            if(fuelRatio_[ne]>0.0 && !stSpecies_[ns]) n0Fuel[ne] += elements_[ne][ns] * X0[ns];
            if(currentSpecie != oxidizer_ && !stSpecies_[ns]) n0Fuel[ne] += elements_[ne][ns] * X0[ns]; 
	    }
    }

    //- In 1 mol of gas (sp), additional O-element required to convert all radicals into its product
    scalar nOxLoss = 0.0;
    forAll(elements_,ne) nOxLoss += n0Fuel[ne] * fuelRatio_[ne];
    nOxLoss /= oxRatio_; //- nOxLoss: number of elements, oxRatio: number of ox elemnts per ox molecule

    // In 1 mol of gas (sp), O2 left after the conversion 
    nOxNet_ = nOx - nOxLoss;
    if(printOut_)  Info << "nOxNet_= " << nOxNet_ << " = " << nOx << " - " << nOxLoss << " (nOx-nOxLoss)" << endl;
}

template<class CloudType>
void Foam::EquilibriumApproximation<CloudType>::solveTable
(
    const scalarField Y1,
    const scalar f0
)
{
    scalar MW1 = this->owner().pope().composition().mixtureMW(Y1);
    scalarField X1 = this->owner().pope().composition().X(Y1);
    scalarField n1;
    n1.setSize(nElements_,0.0);
    YSt_.setSize(this->owner().pope().composition().componentNames().size(),0.0);
//    label nElementsOx = 0;

    //- Populating n1 group
    forAll(X1,ns)
    {
        forAll(elements_,ne)
        {
            n1[ne] += elements_[ne][ns] * X1[ns];
//            if(this->owner().pope().composition().componentNames()[ns] == oxidizer_) nElementsOx += elements_[ne][ns];
        }
    }

    //- To consume 1 mol of fuel, quantity of pure oxidizer required
    scalar nOxReq = 0.0;
    forAll(elements_,ne) nOxReq += n1[ne] * fuelRatio_[ne];
    nOxReq /= oxRatio_;

    //- To consume 1 mol of fuel, quantity of gas (sp) required
    scalar n0Req = nOxReq / nOxNet_;
    if(printOut_) Info << "n0Req=nOxReq/nOxNet: " << n0Req << " = " << nOxReq << "/" << nOxNet_ << endl;

    //- Calculate quantity of products at stoichimetry
    //  Balanced reaction: 
    //  1 {fuel} + n0Req {SP} ==> nProd[CO2] + nProd[H2O] + nProd[N2] (example for stoichometry with CO2, H2O and N2) 
    scalarField nProd;
//    label nOxElementsProd = 0;
//    label nElementsProd = 0;
    nProd.setSize(nSpeciesSt_,0.0);
    scalar mProd = 0.0;

    forAll(nProd,np)
    {
        if(printOut_) Info << speciesSt_[np] << ": ";
        forAll(elements_,ne)
        {
            if(fuelRatio_[ne]>=0.0 && elements_[ne][stThermoId_[np]] > 0)
            {
                nProd[np] += (n1[ne] + n0Req * n0_[ne])/elements_[ne][stThermoId_[np]];
                if(printOut_) Info << " +(" << n1[ne] <<" + " << n0Req << "+" << n0_[ne] << ")/" << elements_[ne][stThermoId_[np]] << ") [" << (n1[ne] + n0Req * n0_[ne])/elements_[ne][stThermoId_[np]] << "] " << elementNames_[ne] << " ";
            }
        }
        //- Compute local stoich. composition (final product)
        YSt_[stThermoId_[np]] = nProd[np] * this->owner().pope().composition().molWt(stThermoId_[np]);
	    mProd += YSt_[stThermoId_[np]];
    }
    YSt_ /= mProd;
    if(printOut_) Info << "\nmProd: " << mProd << endl;

    //- Compute stoich. mixture fraction
    scalar fStTmp = 1.0 * MW1 / (1.0 * MW1 + n0Req * MW0_);
    fSt_ = fStTmp + f0 - fStTmp * f0; //- shift fSt in case that fGas > 0

    if(printOut_) 
    {
        Info << "MW1= " << MW1 << ", n0Req= " << n0Req << ", MW0_= " << MW0_ << ", fSt:  " << fSt_ << endl;
        Info << "YSt: " ;
        forAll(YSt_, nYSt)
        {
            if(YSt_[nYSt]>0.0) Info << YSt_[nYSt] << " " << this->owner().pope().composition().componentNames()[nYSt] << ", " ;
        }
        Info <<endl;
    }
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
    this->setBoundaries(f0,1.0);

    if(envelopeSwitch_ && fSurfOld_ >= fSt_ && fSurfOld_ != 2.0) //- fSurfOld_ != 2.0: fSurf computation in first iteration
    {                                                            //  always based on pure mixing (no envelope)
        if(printOut_) Info << "fSurface CASE I: fSt_= " << fSt_ << ", fSurfOld=" << fSurfOld_;
        if(Y1Liq != Y1LiqOld_) // -in case for differential evaporation YFuelVap can change during iteration,
        {                      //  for non diff evap not -> solve table only once
            this->solveTable(this->Y1(Y1Liq),this->f0_);
            Y1LiqOld_ = Y1Liq;
        }
        fSurf = this->fEnvelope(YSt_,Y1Liq,fSt_,this->f1_,YSurf);
    }
    else
    {
        if(printOut_) Info << "fSurface CASE II: fSt_= " << fSt_ << ", fSurfOld=" << fSurfOld_;
        fSurf = this->fNoEnvelope(Y0,Y1Liq,YSurf);
    }
    if(printOut_ && envelopeSwitch_) Info << "fSurf: " << fSurf << endl;

    return fSurf;
}


template<class CloudType>
Foam::scalar Foam::EquilibriumApproximation<CloudType>::fEnvelope
(
    const scalarField Y0,
    const scalarField Y1,
    const scalar f0,
    const scalar f1,
    const scalarField YX
) const
{
    scalarField linkFG = this->owner().linkFG();
    scalar linkFGSize = linkFG.size();
    scalar fX = -1.0;
    scalar fSum = 0.0;

    forAll(linkFG,nfs)
    {
        fX = (f1 - f0) / (Y1[nfs] - Y0[linkFG[nfs]] + VSMALL) * (YX[nfs] - Y0[linkFG[nfs]]) + f0;
        fSum += fX;
    }
    fX = fSum / linkFGSize;

    return fX;
}

template<class CloudType>
Foam::scalarField Foam::EquilibriumApproximation<CloudType>::YSurface
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar fSurf
)
{
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YSurf;
    YSurf.setSize(Y1.size(),0.0);

    if(envelopeSwitch_ && fSurf >= fSt_ && fSurfOld_ != 2.0) YSurf = this->YEnvelope(YSt_,Y1,fSt_,this->f1_,fSurf);
    else YSurf = this->YNoEnvelope(Y0,Y1,fSurf);

//**//    fSurfOld_ = fSurf[this->owner().linkFG().size()-1];
    fSurfOld_ = fSurf;

    return YSurf;
}

template<class CloudType>
Foam::scalarField Foam::EquilibriumApproximation<CloudType>::YFilm
(
    const scalarField Y0,
    const scalarField Y1Liq,
    const scalar fFilm
) const
{   
    scalarField Y1 = this->Y1(Y1Liq);
    scalarField YFilm;
    YFilm.setSize(Y1.size(),0.0);

    if     (envelopeSwitch_ && fSurfOld_ > fSt_ && fFilm >= fSt_) YFilm = this->YEnvelope(YSt_,Y1,  fSt_,this->f1_,fFilm);
    else if(envelopeSwitch_ && fSurfOld_ > fSt_ && fFilm < fSt_)  YFilm = this->YEnvelope(Y0,  YSt_,this->f0_,fSt_,fFilm);
    else                                                          YFilm = this->YNoEnvelope(Y0,  Y1, fFilm); 

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
    scalar sumY = 0.0;

    forAll(YX, ns)
    {
        YX[ns] = (Y1[ns] - Y0[ns]) / (f1 - f0 + VSMALL) * (fX - f0) + Y0[ns];
        sumY += YX[ns];
    }
    YX /= sumY;

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

    haSurf = this->haNoEnvelope(YSurf,p,TD); //- only for droplet evaporation, needs to be changed for coal

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

    if(printOut_) Info << "fSurf= " << fSurf << ", fSt= " << fSt_ << ", f0= " << this->f0_ << ", f1= " << this->f1_ << endl;
    haX = (ha1 - ha0) / (this->f1_ - this->f0_) * (fSurf - this->f0_) + ha0;

    return haX;
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

    TSurf = this->TNoEnvelope(TD);  //- only for droplet evaporation, needs to be changed for coal

    return TSurf;
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
    scalar sumY = 0.0;
    forAll(YSurf,ns) sumY += YSurf[ns];

    TSurf = this->owner().pope().composition().particleMixture(YSurf).THa(haSurf, p, T0);

    return TSurf;
}

// ************************************************************************* //
