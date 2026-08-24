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

#include "KernelEstimation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
void Foam::KernelEstimation<CloudType>::buildLESParticleList()
{
    const HashTable<label, word>& fieldIndexes = this->XiC().cVarInXi();
    const HashTable<label, word>& pfieldIndexes = this->XiC().cVarInXiC();

    LESList_.setSize(mesh_.cells().size());

    forAll(mesh_.cells(),celli)
    {
        densParticle p(numP_);

        p[nI_]  = celli;
        p[nx_]  = mesh_.C()[celli][0];
        p[ny_]  = mesh_.C()[celli][1];
        p[nz_]  = mesh_.C()[celli][2];
        p[nwt_] = 0.;
        p[nT_]  = 0.;

        //- Radial distance to origin
        p[nDist_] = mag(mesh_.C()[celli]);

        //- Species
        label n = nYs_;

        forAllIters(Yindexes_, Yi)
        {
            p[n] = 0.;//YEqvETarget[*Yi][celli];
            n++;
        }

        forAllConstIter(wordList,this->XiCNames(), nameI)
        {
            //p[nXiCs_ +  pfieldIndexes[*nameI]] =
            //    this->XiC().Vars(fieldIndexes[*nameI]).field()[celli];
	    const label slot = pfieldIndexes[*nameI];

            // In phi-degree mode the carrier slot holds the per-cell projected
            // phi-degree instead of the registered Eulerian coupling value.
            p[nXiCs_ + slot] =
                (phiModEnabled_ && slot == condSlotXiC_)
              ? phiModCell_()[celli]
              : this->XiC().Vars(fieldIndexes[*nameI]).field()[celli];
        }

        LESList_[celli] = p;
    }
}


template <class CloudType>
void Foam::KernelEstimation<CloudType>::buildPhiModCell()
{
    // Project the (flagged) particles' modified progress variable phi-degree
    // onto the mesh as the per-cell conditioning coordinate. phi-degree is a
    // particle-only quantity, so without this projection the kernel would have
    // no cell value to condition on. Super-cells are used (as in
    // ParticleInCell) so the field is defined wherever a flagged particle
    // exists in the super-cell; cells with no flagged particle keep the
    // sentinel -1 and are skipped by the [fLow,fHigh] gate in computeTargets.
    phiModCell_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiModCell",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("phiModCell", dimless, -1.0)
        )
    );

    scalarField& pm = phiModCell_->primitiveFieldRef();

    const bool filterFlagged = this->owner().secondCondMixingEnabled();

    const auto& pManager = this->owner().pManager();

    auto cellParticlesInSuperCell =
        pManager.getParticlesInSuperCellList(this->owner());

    const auto& superCellCells = pManager.cellsInSuperCell();

    for
    (
        label superCelli = 0;
        superCelli < pManager.nSuperCells();
        ++superCelli
    )
    {
        const DynamicList<particleType*>& cellParticles =
            cellParticlesInSuperCell[superCelli];

        scalar sumWt    = 0.0;
        scalar sumPhiWt = 0.0;

        for (particleType* pPtr : cellParticles)
        {
            if (filterFlagged && pPtr->secondCondFlag() != 1) continue;

            const scalar w = pPtr->wt();
            sumWt    += w;
            sumPhiWt += w * pPtr->phiModified();
        }

        if (sumWt > SMALL)
        {
            const scalar mean = sumPhiWt/sumWt;
            for (const label celli : superCellCells[superCelli])
                pm[celli] = mean;
        }
    }
}


template <class CloudType>
void Foam::KernelEstimation<CloudType>::buildParticleList()
{
    particleList_.clear();


    const bool filterFlagged = this->owner().secondCondMixingEnabled();

    // Down-sampling (to ~0.5*nCells particles) is only applied to the full
    // cloud. With the second-conditioning subset active EVERY flagged particle
    // is kept: the subset is already reduced and thinning it would re-introduce
    // the sparse coupling holes we are trying to avoid.
    const scalar nParticles = this->owner().size();
    const scalar nLESCells  = mesh_.cells().size();
    const bool   subSample  = (!filterFlagged) && (nParticles > nLESCells);
    const scalar prob       = subSample ? 1. - (0.5*nLESCells/nParticles) : 0.0;


    forAllIters(this->owner(), iter)
    {
		// When second conditioning is active, exclude non-subset particles —
        // their Y and T are stale and must not bias the Eulerian target fields.
        //if (this->owner().secondCondMixingEnabled() && iter().secondCondFlag() != 1)
        if (filterFlagged && iter().secondCondFlag() != 1)
            continue;
		
        densParticle p(numP_);

        p[nI_] = 0.;

        p[nx_] = iter().position().x();

        p[ny_] = iter().position().y();

        p[nz_] = iter().position().z();

        p[nwt_] = iter().m();

        p[nT_] = iter().T();

        p[nDist_] = VGREAT;

        //- Species
        label n = 0;
        forAllIter(labelList,Yindexes_, Yi)
        {
            p[nYs_ + n] = iter().Y()[*Yi];

            //- Find min and max of particle for each YEqvE
            maxVal_[n] = max(maxVal_[n],iter().Y()[*Yi]);
            minVal_[n] = min(minVal_[n],iter().Y()[*Yi]);
            n++;
        }

        //- Coupling(state) Variables
        forAll(iter().XiC(), j)
        {
            //p[nXiCs_ + j] = iter().XiC()[j];
            //maxXiCVal_[j] = max(maxXiCVal_[j],iter().XiC()[j]);
            //minXiCVal_[j] = min(minXiCVal_[j],iter().XiC()[j]);
	    // In phi-degree mode the carrier slot holds phiModified() so the
            // kd-tree conditions on the modified progress variable.
            const scalar val =
                (phiModEnabled_ && j == condSlotXiC_)
              ? iter().phiModified()
              : iter().XiC()[j];

            p[nXiCs_ + j] = val;
            maxXiCVal_[j] = max(maxXiCVal_[j], val);
            minXiCVal_[j] = min(minXiCVal_[j], val);
        }

        //- Min and Max of temperature
        maxVal_[numYEqv_] = max(maxVal_[numYEqv_],iter().T());
        minVal_[numYEqv_] = min(minVal_[numYEqv_],iter().T());

        //if(nParticles > nLESCells)
	if (subSample)
        {
            if (this->owner().rndGen().Random() >= prob)
                particleList_.append(p);
        }
        else
            particleList_.append(p);
    }
}

template <class CloudType>
void Foam::KernelEstimation<CloudType>::computeTargets
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{
    label N2Index = 0;

    volScalarField Yt(0.0*YEqvETarget[0]);

    //- Indicator (used to make zero source terms)
    this->Indicator() = Yt;

    TEqvETarget = this->owner().T();

    scalarField dTEqvETarget(TEqvETarget.size(),0.0);

    List<scalarField> dYEqvETarget;

    forAll(YEqvETarget, specieI)
    {
        if (!this->solveEqvSpecie()[specieI])
            YEqvETarget[specieI].primitiveFieldRef() = 0.0;
        else
        {
            YEqvETarget[specieI] = this->owner().composition().carrier().Y()[specieI];
           dYEqvETarget.append(scalarField(TEqvETarget.size(),0.0));

            Yt += YEqvETarget[specieI];
        }

        if (YEqvETarget[specieI].name()==inertSpecie_)
            N2Index = specieI;
    }

    //- Index of coupling variable used
        //const label CVIndexinXiC = nXiCs_ + XiC.cVarInXiC()[this->cVarName()];
    //- Index of conditioning variable used (carrier slot; in phi-degree mode
    //- the first coupling variable's slot carries the conditioning value)
    const label CVIndexinXiC = nXiCs_ + condSlotXiC_;
    const label CVIndexinXi  = condSlotXi_;

    //- Min and Max values of the coupling variable
    //const label CVIndexinXi = XiC.cVarInXi()[this->cVarName()];
    //- LES field of the conditioning variable: the projected phi-degree field
    //- in phi-degree mode, otherwise the registered Eulerian coupling field.
    const volScalarField& f =
        phiModEnabled_ ? phiModCell_() : XiC.Vars(CVIndexinXi).field();

    //scalar minf = min(XiC.Vars(CVIndexinXi).field().primitiveField());
    //scalar maxf = max(XiC.Vars(CVIndexinXi).field().primitiveField());
    //- Min and Max values of the conditioning variable
    scalar minf = min(f.primitiveField());
    scalar maxf = max(f.primitiveField());

    //scalar minz = minXiCVal_[XiC.cVarInXiC()[this->cVarName()]];
    //scalar maxz = maxXiCVal_[XiC.cVarInXiC()[this->cVarName()]];
    scalar minz = minXiCVal_[condSlotXiC_];
    scalar maxz = maxXiCVal_[condSlotXiC_];

    if (debug_)
        Info << "max and min values of couplingVar in LES are: "
             << maxf << " and "<< minf
             << "\nmax and min values of couplingVar in particles are: "
             << maxz << " and "<< minz
             << "\nmax and min values of YEqv in particles are: \n"
             << maxVal_ << nl << minVal_ << endl;

    //- Pointer Lists
        //- LES Particles List

        List<densParticle*> LESPtrList(LESList_.size());

        forAll(LESList_,ii)
               LESPtrList[ii] = &LESList_[ii];

    //- Iterators for sorting
    List<densParticle*>::iterator iterL2 = LESPtrList.begin();
    List<densParticle*>::iterator iterU2 = LESPtrList.end();

    if (debug_)
        Info << "Iterators initialized" << endl;

    //- Sort particleList

        labelList dims(4);
        dims[0] = nx_;
        dims[1] = ny_;
        dims[2] = nz_;
        dims[3] = CVIndexinXiC; 

        scalarList wts(4);
        wts[0] = 1.;
        wts[1] = 1.;
        wts[2] = 1.;
        wts[3] = fm_;

        // Convert the particleList to a List<List<scalar>> with only the 
        // dimensions of dims for kdTree 
        // create the list
        List<List<scalar>> pList(particleList_.size());

        forAll(pList,i)
        {
            List<scalar> temp(4);
            for (int j=0; j < 4; ++j)
                temp[j] = particleList_[i][dims[j]];
            pList[i] =std::move(temp);
        }

        //- Construct kd-Tree of particle list
        //kdTree<List<scalar>> particleTree(pList,wts); // 4 dimensions for kdtree x,y,z,XiC
	const bool haveParticles = !pList.empty();

        autoPtr<kdTree<List<scalar>>> particleTreePtr;
        if (haveParticles)
            particleTreePtr.reset
            (
                new kdTree<List<scalar>>(pList, wts) // 4 dims: x,y,z,condVar
            );

        //Sort LES list based on distance from origin
        std::sort(iterL2,iterU2,lessArg(nDist_));

    if (debug_)
        Info << "Lists sorted" << endl;

    //- Field with cell centres
    vectorField cellCentres = mesh_.C().internalField();

    //- LES field of coupling variable
    //const volScalarField& f = XiC.Vars(CVIndexinXi).field();//

    volVectorField gradf = fvc::grad(f);

    if (debug_)
        Info << "Start to compute cell averages" << endl;

    label computedLES = 0;// counter for kernel computed cells

    //- Index of last cell computed with kernel
    label lastCCell = (*LESPtrList.first())[nI_];

    scalar DELTAd = 0.0;
    scalar DELTAf = 0.0;

    //- Start computation of Cell target values
    forAllIter(DynamicList<densParticle*>,LESPtrList,LESi)
    {
        if (!haveParticles)
            break;

	scalar fLES  = (**LESi)[CVIndexinXiC];

        label celli  = (**LESi)[nI_];

        scalar f0LES = f[lastCCell];
        scalar df    = fLES - f0LES;

        //- Do not compute cell value, source term is zero
        //- Used to avoid biased kernel at boundaries in coupling space
        //- It also helps to speed up computations
        if ((fLES <= fLow_) || (fLES >= fHigh_))
            continue;

        //- Cell centre
        vector cCentre  = cellCentres[celli];
        vector c0Centre = cellCentres[lastCCell];

        //- Physical distance
        scalar dd = mag(cCentre - c0Centre);

        //- compute using kernel if cell is far in coupling space (df>= DELTAf)
        //- or physical space (dd>= DELTAd)
        if (mag(df) >= DELTAf || dd >= DELTAd || computedLES == 0)
        {
            //- Find the k-nn particles to compute kernel !!!!!!!!
            //label nn = 20;//50;
	    const label nn = nNearest_;
            scalarList qv(4); //coupling var + 3 dimensions

            qv[0] = cCentre[0];
            qv[1] = cCentre[1];
            qv[2] = cCentre[2];
            qv[3] = fLES;

            //auto result = particleTree.nNearest
	    auto result = particleTreePtr().nNearest
            (
                qv,             // query vector (x,y,z,XiC)
                nn              // number of nearest neighbours
            );

	    // The tree returns at most nn neighbours (fewer when the flagged
            // subset is smaller than nn); guard against an empty/short result
            // before sizing the kernel from result.begin() and looping.
            if (result.empty())
                continue;

            const label nFound = min(nn, label(result.size()));

            //- Define kernel radius in physical space based on k-NN

            vector disVector;
            disVector[0] = result.begin()->disReal[0];
            disVector[1] = result.begin()->disReal[1];
            disVector[2] = result.begin()->disReal[2];

            scalar rMax = min(max(0.0,mag(disVector)),rMaxMax_);

            //- Start to compute target values
            scalar sumWt   = 0.0;
            scalar sumdWt  = 0.0;
            scalar sumWtT  = 0.0;
            scalar sumdWtT = 0.0;
            scalarField sumWtYi(Yindexes_.size(),0.0);
            scalarField sumdWtYi(Yindexes_.size(),0.0);

            //- Loop over k-NN Paricles

            scalar h  = 0.25*fm_;
            scalar h2 = 0.5*(rMax+SMALL);

            scalar alpha   = 1.0/6.0/h;
            scalar alpha2  = alpha/(0.5*h);

            //for(label nni = 0; nni < nn; nni++)
	    for(label nni = 0; nni < nFound; nni++)
            {
                const densParticle& p = particleList_(result[nni].idx);

                //- Distance in mix fraction
                scalar cpfDistance = -result[nni].disReal[3];

                scalar NfDistance = cpfDistance/h; //Normalized distance

                vector disXP;
                disXP[0] = result[nni].disReal[0];
                disXP[1] = result[nni].disReal[1];
                disXP[2] = result[nni].disReal[2];

                scalar NDistance  = mag(disXP)/h2; //////!!!!!!!!

                //- 1-D cubic smoothing kernel weight
                //- J J Monaghan Rep. Prog. Phys.68 (2005) 1703–1759
                //- doi:10.1088/0034-4885/68/8/R01

                scalar  IDWf   = 0.0;
                scalar dIDWf   = 0.0;

                scalar q    = mag(NfDistance);
                scalar sign = (NfDistance > 0.0) ? 1.0 : ((NfDistance < 0.0) ? -1.0 : 0.0);

                if (q >= 0.0 && q < 1.0)
                {
                    IDWf = alpha*(pow(2.0 - q, 3.0) - 4.0 * pow(1.0 - q, 3.0));
                   dIDWf = -sign*alpha2*(3.0*sqr(2.0 - q) - 12.0*sqr(1.0 - q));
                }
                else if (q >= 1.0 && q < 2.0)
                {
                    IDWf = alpha*(pow(2.0 - q, 3.0));
                   dIDWf = -sign*alpha2*(3.0*sqr(2.0 - q));
                }

                //- 3-D cubic smoothing kernel weight
                //- J J Monaghan Rep. Prog. Phys.68 (2005) 1703–1759
                //- doi:10.1088/0034-4885/68/8/R01

                scalar IDWd   = 0.0;
                scalar alphad = 1.0/(4.0*pi*pow(h2,3.0));

                if (NDistance < 1.0)
                    IDWd = alphad*(pow(2.0 - NDistance,3.0) - 4.0 * pow(1.0 - NDistance,3.0));
                else if ((NDistance >= 1.0) && (NDistance < 2.0))
                    IDWd = alphad*(pow(2.0 - NDistance,3.0));

                scalar mWti = IDWd * p[nwt_];
                scalar  Wti = mWti *  IDWf;
                scalar dWti = mWti * dIDWf;

                sumWt   = sumWt   +  Wti;
                sumdWt  = sumdWt  + dWti;

                sumWtT  = sumWtT  + (p[nT_]  *  Wti);
                sumdWtT = sumdWtT + (p[nT_]  * dWti);

                label ii = 0;
                for (label n = nYs_; n <= nYe_ ; n++)
                {
                    sumWtYi[ii]  = sumWtYi[ii]  + (p[n] *  Wti);
                    sumdWtYi[ii] = sumdWtYi[ii] + (p[n] * dWti);
                    ii++;
                }
            }

            if (sumWt < VSMALL)
                continue;

            Yt[celli]      = 0.0;
            label I        = 0;
            scalar maxGrad = 0.0;

            forAllIter(labelList, Yindexes_, specieI)
            {
                scalar Yi = sumWtYi[I]/sumWt;

                YEqvETarget[*specieI][celli] = Yi;
               dYEqvETarget[I][celli] = (sumdWtYi[I] - (Yi * sumdWt))/sumWt;

                //- Find maximum gradient to determine distance allowed to extrapolate
                scalar diffMag  = mag(maxVal_[I]-minVal_[I]);
                scalar Ngrad    = (diffMag > ROOTVSMALL) ? diffMag/h : 1e30;
                scalar maxGradi = mag(dYEqvETarget[I][celli]/Ngrad);

                if (maxGradi>maxGrad)
                    maxGrad = maxGradi;

                if (*specieI != N2Index)
                    Yt[celli] = Yt[celli] + YEqvETarget[*specieI][celli];

                I++;
            }

            TEqvETarget[celli]  = sumWtT/sumWt;
            dTEqvETarget[celli] = (sumdWtT - (TEqvETarget[celli] * sumdWt))/sumWt;

            //- Find maximum gradient to determine distance allowed to extrapolate
            scalar diffMag  = mag(maxVal_[I]-minVal_[I]);
            scalar Ngrad    = (diffMag > ROOTVSMALL) ? diffMag/h : 1e30;
            scalar maxGradi = mag(dTEqvETarget[celli]/Ngrad);

            if (maxGradi>maxGrad)
            maxGrad = maxGradi;

            DELTAf = C2_*h/(maxGrad+ROOTVSMALL);// define C2 !!! Needs checking

            if (maxGrad >= 10.0)
                DELTAf = 0.0;
            else
                DELTAf = min(dfMax_,DELTAf);// define DELTA_MAX

            scalar gradFi = cmptMax(cmptMag(gradf[celli]));

            DELTAd = DELTAf/(gradFi +ROOTVSMALL);

            this->Indicator()[celli] = 1.0;

            if (debug_ && (TEqvETarget[celli]>3000.0 || TEqvETarget[celli] < 290.0))
                Pout << "Target Temp is: " << TEqvETarget[celli] << endl;

            lastCCell = celli;

            computedLES++;
        }
        else
        {
        //- Compute target value using Linear Extrapolation from
        //- last cell computed
            Yt[celli] = 0.0;

            label I = 0;
            forAllIter(labelList,Yindexes_,specieI)
            {
                scalar Yi = YEqvETarget[*specieI][lastCCell];
                scalar deltaYi = dYEqvETarget[I][lastCCell]*df;

                YEqvETarget[*specieI][celli] = Yi + deltaYi;

                if (*specieI != N2Index)
                    Yt[celli] = Yt[celli] + YEqvETarget[*specieI][celli];

                I++;
            }

            scalar deltaT = dTEqvETarget[lastCCell]*df;

            TEqvETarget[celli] = TEqvETarget[lastCCell] + deltaT;

            this->Indicator()[celli] = 1.0;

            if (debug_ && (TEqvETarget[celli]>2200.0 || TEqvETarget[celli] < 290.0))
            {
                    Info << "fLES: " << fLES << endl;
                    Info << "Target Temp is: " << TEqvETarget[celli] << endl;
            }
        }
    }

    //- Make species add up to one
    YEqvETarget[N2Index] = scalar(1) - Yt;
    YEqvETarget[N2Index].max(0.0);

    // Coupling-coverage diagnostic: fraction of cells that received a Lagrangian
    // target (Indicator==1). With the kernel estimator this should stay high even
    // for the sparse flagged subset, unlike ParticleInCell which leaves a hole
    // wherever a flagged particle is absent from the cell.
    {
        const scalar nCov =
            returnReduce(sum(this->Indicator().primitiveField()), sumOp<scalar>());
        const label nTot =
            returnReduce(this->Indicator().size(), sumOp<label>());
        const label nP = returnReduce(particleList_.size(), sumOp<label>());

        Info<< "KernelEstimation coupling: " << label(nCov) << "/" << nTot
            << " cells covered (" << 100.0*nCov/max(nTot, 1) << "%), "
            << nP << " particles in kernel list" << endl;
    }

    if (debug_)
    {
        Pstream::gather(computedLES,sumOp<label>());
        Info << "Target Values Computed \nThe number of cells with a computed "
             << "target value using kernel estimator is: " << computedLES << endl;
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KernelEstimation<CloudType>::KernelEstimation
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    ThermoPhysicalCouplingModel<CloudType>(dict,owner,typeName,Xi),

    mesh_(owner.mesh()),

    Yindexes_(),

    numYEqv_(0),

    nI_(0),

    nx_(nI_ + 1),

    ny_(nx_ + 1),

    nz_(ny_ + 1),

    nwt_(nz_ + 1),

    nT_(nwt_ + 1),

    nDist_(nT_ + 1),

    nYs_(nDist_ + 1),

    nYe_(nYs_ + 1),

    nXiCs_(nYe_ + 1),

    nXiCe_(nXiCs_ + 1),

    numP_(0),

    fm_(readScalar(this->coeffDict().lookup("fm"))),

    rMaxMax_(this->coeffDict().lookupOrDefault("rMax", 1.0e9)),

    nNearest_(this->coeffDict().lookupOrDefault("nNearest", label(20))),

    particleList_(),

    LESList_(),

    maxVal_(),

    minVal_(),

    maxXiCVal_(),

    minXiCVal_(),

    fLow_(readScalar(this->coeffDict().lookup("fLow"))),

    fHigh_(readScalar(this->coeffDict().lookup("fHigh"))),

    inertSpecie_(this->owner().thermo().lookup("inertSpecie")),

    dfMax_(readScalar(this->coeffDict().lookup("dfMax"))),

    C2_(readScalar(this->coeffDict().lookup("C2"))),

    //debug_(this->coeffDict().lookupOrDefault("debug",false)) //Define from dictionary
    debug_(this->coeffDict().lookupOrDefault("debug",false)), //Define from dictionary

    phiModEnabled_(false),

    condSlotXiC_(0),

    condSlotXi_(0),

    phiModCell_(nullptr)
{
    forAll(this->solveEqvSpecie(), i)
    {
        if (this->solveEqvSpecie()[i])
            Yindexes_.append(i);
    }

    numYEqv_ = Yindexes_.size();
    numXiC_  = this->XiCNames().size();

    nYe_    = nYs_ + numYEqv_ - 1;
    nXiCs_  = nYe_ + 1;
    nXiCe_  = nXiCs_ + numXiC_ - 1;
    numP_   = nXiCe_ + 1;

    maxVal_.setSize(numYEqv_+1,0.0);
    minVal_.setSize(numYEqv_+1,GREAT);

    maxXiCVal_.setSize(numXiC_,0.0);
    minXiCVal_.setSize(numXiC_,GREAT);

    if (debug_)
        Info << "The number of equivalent species is: " << numYEqv_
             << " and the number of properties is: " << numP_ << nl
             << "The indices of the equivalents species are: "
             << Yindexes_ << endl;
    Info << "maximum rMax is: " << rMaxMax_ << endl;

    // Optional conditioning on the modified progress variable phi-degree.
    // phi-degree (phiModified) is a particle-only quantity with no Eulerian
    // field, so the first registered coupling variable's array slot is reused
    // as the carrier for the conditioning value and a per-cell phi-degree is
    // projected from the flagged particles each step (buildPhiModCell()).
    phiModEnabled_ = (this->cVarName() == "phiModified");

    if (phiModEnabled_ && this->XiCNames().empty())
        FatalErrorInFunction
            << "condVariable 'phiModified' requires at least one registered "
            << "coupling variable to carry the conditioning slot." << nl
            << exit(FatalError);

    const word condName =
        phiModEnabled_ ? this->XiCNames()[0] : this->cVarName();

    condSlotXiC_ = this->XiC().cVarInXiC()[condName];
    condSlotXi_  = this->XiC().cVarInXi()[condName];

    if (phiModEnabled_)
    {
        Info<< "KernelEstimation: conditioning on phiModified (phi-degree), "
            << "carrier coupling slot '" << condName << "'" << endl;

        // phi_deg = phi*exp(beta*omega) is unbounded above even though phi is
        // clamped to [0,1]: omega ~ N(0,1), so exp(beta*omega) exceeds 1 for
        // half the particles. fLow/fHigh are normally calibrated for a mixture
        // fraction in [0,1], and computeTargets() skips every cell whose
        // conditioning value falls outside them -- leaving Indicator = 0 and
        // switching the relaxation source off there. Warn rather than silently
        // adjusting fHigh, and point at the coverage diagnostic.
        if (fHigh_ <= 1.0)
        {
            WarningInFunction
                << "condVariable is 'phiModified' but fHigh = " << fHigh_
                << " is calibrated for a variable bounded by 1." << nl
                << "phiModified = phi*exp(beta*omegaOU) is NOT bounded above, "
                << "so cells above fHigh will be dropped from the coupling "
                << "(Indicator = 0)." << nl
                << "Check the 'KernelEstimation coupling: N/M cells covered' "
                << "line and the phiModified range reported by "
                << "secondCondMMCcurl, and raise fHigh accordingly."
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KernelEstimation<CloudType>::~KernelEstimation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KernelEstimation<CloudType>::EqvETargetValues
(
    const mmcVarSet& XiC,
    PtrList<volScalarField>& YEqvETarget,
    volScalarField& TEqvETarget
)
{
    //- In phi-degree mode, project the flagged particles' phi-degree onto the
    //- mesh first: it is the per-cell conditioning coordinate the kernel needs.
    if (phiModEnabled_)
        buildPhiModCell();

    //- Create lists
    buildParticleList();

    buildLESParticleList();

    //- Gather particles of neighbour processors to avoid bias at 
    //- processor boundaries.

    labelList pNeighbours;

    forAll(mesh_.boundaryMesh(),patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& ppp = 
                refCast<const processorPolyPatch>(pp);

            pNeighbours.append(ppp.neighbProcNo());
        }
    }

    //- Use particles in processor plus neighbouring processors
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (label procI : pNeighbours)
        {
            UOPstream toBuffer(procI, pBufs);
            toBuffer << particleList_;
        }
        
        pBufs.finishedSends();

        for (label procI : pNeighbours)
        {
            DynamicList<densParticle> tempList;
            UIPstream fromBuffer(procI, pBufs);
            fromBuffer >> tempList;
            // Add particles of processor to current list
            particleList_.append(tempList);
        }
    }

    computeTargets(XiC,YEqvETarget,TEqvETarget);

    LESList_.clearStorage();
    particleList_.clearStorage();
}

// ************************************************************************* //


