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

#include "ransMMCcurl.H"
#include "fvMesh.H"
#include "StochasticLib.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template <class CloudType>
List<scalar> Foam::ransMMCcurl<CloudType>::getXiNormalisation()
{
    // dictionary to read the normalisation parameters for 
    // the reference variables 
    const dictionary XiDict(this->coeffDict().subDict("Xim_i"));

    Info << nl << "The Xim_i parameters are: "<< XiDict << endl;

    List<scalar> Xii(this->numXiR());

    label i=0;
    for (const word& refVarName :this->XiRNames())
    {
        Xii[i++] = readScalar(XiDict.lookup(refVarName+"m"));
    }

    return Xii;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ransMMCcurl<CloudType>::ransMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    ransMixParticleModel<CloudType>(dict,owner, typeName, Xi),

    C_phi_(readScalar(this->coeffDict().lookup("C_phi"))),

    rt_(readScalar(this->coeffDict().lookup("r_t"))),

    bo_(readScalar(this->coeffDict().lookup("b_o"))),

    Xi_i_(getXiNormalisation()),

    XiR_(Xi)
{}


template <class CloudType>
Foam::ransMMCcurl<CloudType>::ransMMCcurl
(
    const ransMMCcurl<CloudType>& cm
)
:
    ransMixParticleModel<CloudType>(cm),

    C_phi_(readScalar(this->coeffDict().lookup("C_phi"))),

    rt_(readScalar(this->coeffDict().lookup("r_t"))),

    bo_(readScalar(this->coeffDict().lookup("b_o"))),

    Xi_i_(getXiNormalisation()),

    XiR_(cm.XiR_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ransMMCcurl<CloudType>::~ransMMCcurl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ransMMCcurl<CloudType>::MixXi
(
    particlePtrList<particleType>& particlePtrList,
    DynamicList<eulerianFieldData>& eulFields
)
{

    scalar deltaT = this->owner().mesh().time().deltaT().value();

    List<scalar> XiSum(Xi_i_.size(),0.0);
    List<scalar> XiErrSum(Xi_i_.size(),0.0);
    scalar wtSum = 0;


    // Calculating the weighted-ensemble-average of v in each PDF cell
    for (label ii=0; ii<particlePtrList.size();ii++)
    {
        particleType& p = particlePtrList[ii];
        
        forAll(Xi_i_,i)
        {
            XiSum[i] = XiSum[i] + p.XiR()[i] * p.wt();
        }
        wtSum = wtSum + p.wt();
    }

    List<scalar> XiAv = XiSum/wtSum;


    // Calculating the weighted-ensemble-average of v in each PDF cell
    for (label ii=0; ii<particlePtrList.size();ii++)
    {

        particleType& p = particlePtrList[ii];

        forAll(Xi_i_,i)
        {
            XiErrSum[i] = XiErrSum[i] + pow((p.XiR()[i] - XiAv[i]),2) * p.wt();
        }
    }

    List<scalar> XiVar = XiErrSum/wtSum;


    for (label ii=0; ii<particlePtrList.size();ii++)
    {
        particleType& p = particlePtrList[ii];
        auto eulerianField = eulFields[ii];
        
        
        scalar tauMix;
        bool mixOn = true;

        if(1.0/eulerianField.tauTurb() < VSMALL)
        {
            tauMix = 1e30;
            mixOn = false;
        }
        else
        {
            tauMix = eulerianField.tauTurb()/C_phi_;
        }


        if (mixOn)
        {
            scalar mixExtent = 1 - exp(-deltaT / (tauMix + VSMALL));

            forAll(Xi_i_,i)
            {
                p.XiR()[i] = p.XiR()[i]  + mixExtent * (XiAv[i] - p.XiR()[i] )
                    +( 
                        this->owner().rndGen().Normal(0,1)
                       * bo_*sqrt(2*XiVar[i]*deltaT/tauMix)
                     );
            }
        }
    }
}


template <class CloudType>
void Foam::ransMMCcurl<CloudType>::mixpair
(
    particleType& p, 
    const eulerianFieldData& pEulField,
    particleType& q,
    const eulerianFieldData& qEulField,
    scalar& deltaT
)
{
    //- If combined weights of p and q > 0
    if(p.wt() + q.wt() > 0)
    {
        scalar tauP;
        scalar tauQ;

        bool mixOn = true;

        if(1.0/pEulField.tauTurb() < VSMALL)
        {
            tauP = 1e30;
            mixOn = false;
        }
        else
        {
            tauP = 2.0*(1.0-pow(rt_,2))*pEulField.tauTurb()/ C_phi_;
        }

        if(1.0/qEulField.tauTurb() < VSMALL)
        {
            tauQ = 1e30;
            mixOn = false;
        }
        else
        {
            tauQ = 2.0*(1.0-pow(rt_,2))*qEulField.tauTurb()/ C_phi_;
        }


        if (mixOn)
        {
            scalar tauMix = max(tauP,tauQ);

            scalar mixExtent = 1.0 - exp(-deltaT / (tauMix + VSMALL));

            // Function cascade to mix scalars
            mixScalarProperties(p,q,mixExtent);
        }
    }
}


template<class CloudType>
void Foam::ransMMCcurl<CloudType>::SmixList
(
    particlePtrList<particleType>& particleList,
    DynamicList<eulerianFieldData>& eulFields
)
{
    // Mix reference variable using new MMC model equation
    if(particleList.size()!=0)
    {
        MixXi(particleList,eulFields);
    }
    
    // Mixing old model
    
    //// Keeping track of indices for kdTreeLikeSearch
    //std::vector<label> L;
    //std::vector<label> U;
    //L.reserve(particleList.size());
    //U.reserve(particleList.size());
    //Info << "kdSearch started"<<endl;
    //kdTreeLikeSearch
    //(
        //1,eulFields.size(),eulFields,L,U
    //);
    //Info << "kdSearch finished"<<endl;
    //for (size_t i=0; i<L.size(); i++)
    //{
        //label p = L[i] - 1;
        
        //label pInd = eulFields[p].particleIndex();

        //label q = L[i];

        //label qInd = eulFields[q].particleIndex();

        //if(U[i] - L[i] < 2)
        //{
            //mixpair
            //(
                //particleList[pInd],eulFields[p],
                //particleList[qInd],eulFields[q],
                //deltaT
            //);
        //}
        //else if(U[i] - L[i] == 2)
        //{
            //label r = L[i] + 1;

            //label rInd = eulFields[r].particleIndex();

            //mixpair
            //(
                //particleList[pInd],eulFields[p],
                //particleList[qInd],eulFields[q],
                //deltaT
            //);

            //mixpair
            //(
                //particleList[qInd],eulFields[q],
                //particleList[rInd],eulFields[r],
                //deltaT
            //);
        //}
    //}
}


template <class CloudType>
void Foam::ransMMCcurl<CloudType>::kdTreeLikeSearch
(
    label l,
    label u,
    DynamicList<eulerianFieldData>& particleList,
    std::vector<label>& L,
    std::vector<label>& U
)
{
//- Break the division if the particle list has length less than 2
    if (u - l <= 2)
    {
        //- Divide particles into groups of two or three
        L.push_back(l);

        U.push_back(u);

        return ;
    }

    label m = (l + u)/2;
    if ( (u - m) % 2 != 0 ) m++;

    auto iterL = particleList.begin();

    auto iterM = particleList.begin();

    auto iterU = particleList.begin();

    std::advance(iterL,l-1);

    std::advance(iterM,m-1);

    std::advance(iterU,u);

    List<scalar> maxInXiR(Xi_i_.size(),0.0);
    List<scalar> minInXiR(Xi_i_.size(),0.0);

    forAll(Xi_i_,i)
    {
        maxInXiR[i] = (*std::max_element(iterL,iterU,lessArg(3 +i))).XiR()[i];
        minInXiR[i] = (*std::min_element(iterL,iterU,lessArg(3 +i))).XiR()[i];
    }

    //- Scaled/stretched distances between Max and Min in each direction
    //- Default is random mixing, overwritten if mixing distances greater than ri or fm
    label ncond = 0;
    scalar disMax = 0;

    forAll(Xi_i_,i)
    {
        scalar disXiR = (maxInXiR[i] - minInXiR[i])/Xi_i_[i];
        if(disXiR > disMax)
        {
            disMax = disXiR;
            ncond = 3+i;
        }
    }

    std::sort(iterL,iterU,lessArg(ncond));

    //- Recursive function calls for lower and upper branches of the particle list
    kdTreeLikeSearch(l,m,particleList,L,U);

    kdTreeLikeSearch(m+1,u,particleList,L,U);
};


template<class CloudType>
const Foam::scalarField Foam::ransMMCcurl<CloudType>::XiR0
(
    label patch,
    label patchFace
)
{

    label nVar = this->numXiR_;
    const mmcVarSet& setOfXi(this->XiR_);

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(nVar,0.0);

    for (const word& varName : this->XiRNames_)
    {
        // Initialisation of the the reference variable from FV fields
        XXi[XiRIndexes[varName]] = setOfXi.Vars(XiIndexes[varName]).field().boundaryField()[patch][patchFace];
    }
    return XXi;

}

template<class CloudType>
const Foam::scalarField Foam::ransMMCcurl<CloudType>::XiR0(label celli)
{
    label nVar = this->numXiR_;
    const mmcVarSet& setOfXi(this->XiR_);

    const labelHashTable& XiIndexes  = setOfXi.rVarInXi();
    const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();

    scalarField XXi(nVar,0.0);

    for (const word& varName : this->XiRNames_)
    {
        // Initialisation of the the reference variable from FV fields
        XXi[XiRIndexes[varName]]=setOfXi.Vars(XiIndexes[varName]).field()[celli];
    }
    return XXi;
}

// NEW FUNCTION TO INITIATE SCALAR PROPERTIES MIXING CASCADE
template<class CloudType>
void Foam::ransMMCcurl<CloudType>::mixScalarProperties
(
    particleType& p,
    particleType& q,
    scalar mixExtent
)
{
    typedef typename CloudType::particleType particleType;

    particleType::mixProperties(p,q,mixExtent);
}


template<class CloudType>
const Foam::mmcVarSet& Foam::ransMMCcurl<CloudType>::XiR() const
{
    return XiR_;
}



// ************************************************************************* //
