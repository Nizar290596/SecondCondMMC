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

#include "nearestNeighbor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, class mmcCloudType>
Foam::nearestNeighbor<CloudType, mmcCloudType>::nearestNeighbor
(
    const dictionary& dict,
    CloudType& owner
)
:
    DropletToMMCModel<CloudType,mmcCloudType>(dict, owner, typeName),
    lambda_x(this->coeffDict().template get<bool>("lambda_x")),
    lambda_f(this->coeffDict().template get<bool>("lambda_f")),
    lambda_T(this->coeffDict().template get<bool>("lambda_T")),

    ri_(this->coeffDict().template get<scalar>("ri")),
    fm_(this->coeffDict().template get<scalar>("fm")),
    Tm_(this->coeffDict().template get<scalar>("Tm")),
    useSPProps_(this->coeffDict().template get<bool>("useSPProps"))
{}


template<class CloudType, class mmcCloudType>
Foam::nearestNeighbor<CloudType, mmcCloudType>::nearestNeighbor
(
    const nearestNeighbor<CloudType, mmcCloudType>& copyModel
)
:
    DropletToMMCModel<CloudType,mmcCloudType>(copyModel.owner_),
    lambda_x(copyModel.lambda_x),
    lambda_f(copyModel.lambda_f),
    lambda_T(copyModel.lambda_T),

    ri_(copyModel.ri_),
    fm_(copyModel.fm_),
    Tm_(copyModel.Tm_),
    useSPProps_(copyModel.useSPProps_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

template<class CloudType, class mmcCloudType>
void Foam::nearestNeighbor<CloudType, mmcCloudType>::buildParticleList
(
    mmcCloudType& mmcCloud
) const
{
    // PreAllocate the 
    particleList_.clear();
    particleList_.reserve(mmcCloud.size());

    forAllIters(mmcCloud,iter)
    {
        vector pos = iter->position();
        scalar fOrZ  = 0;
        if (useSPProps_)
            fOrZ = iter->XiC("z");
        else
            fOrZ = iter->XiR("f");
        particleList_.append
        (
            container
            (
                lambda_x, lambda_f, lambda_T,
                pos.x(),
                pos.y(),
                pos.z(),
                fOrZ,
                iter->T(),
                iter.get()
            )
        );
    }
}


template<class CloudType, class mmcCloudType>
void Foam::nearestNeighbor<CloudType, mmcCloudType>::genKdTree
(
    mmcCloudType& mmcCloud
) const
{
    // Check if the k-d tree has to be re-generated
    label currTimeIndex = this->owner().time().timeIndex();
    if (currTimeIndex != timeIndex_)
    {
        timeIndex_ = currTimeIndex;

        buildParticleList(mmcCloud);

        label size = label(lambda_x)*3+label(lambda_f)+label(lambda_T);

        scalarList wts(size);
        label ind = 0;
        if (lambda_x)
        {
            wts[ind++] = ri_;
            wts[ind++] = ri_;
            wts[ind++] = ri_;
        }
        if (lambda_f)
            wts[ind++] = fm_;
        if (lambda_T)
            wts[ind++] = Tm_;

        kdTree_.reset
        (
            new kdTree<container>
            (
                particleList_,
                wts,
                false,
                false
            )
        );
    }
}


template<class CloudType, class mmcCloudType>
typename Foam::nearestNeighbor<CloudType, mmcCloudType>::mmcParticleType*
Foam::nearestNeighbor<CloudType, mmcCloudType>::locateClosestStochasticParticle
(
    const parcelType& p,    // Particle
    const typename parcelType::trackingData& td,
    mmcCloudType& mmcCloud
) const
{   
    // Get the reference field interpolated from the LES to the particle position
    tetIndices tetIs = p.currentTetIndices();
    const scalar fRef = td.XiRInterp("f").interpolate(p.coordinates(), tetIs);
    const scalar TRef = td.Tc();

    // Generate the k-d tree if necessary
    genKdTree(mmcCloud);

    // Create the reference particle
    // The design of the kdTree class expects here a container object as well
    vector pPos = p.position();
    container convertedP
    (
        lambda_x, lambda_f, lambda_T,
        pPos.x(),
        pPos.y(),
        pPos.z(),
        fRef,
        TRef,
        nullptr
    );


    // Find the closest stochastic particle
    auto res = kdTree_->nNearest(convertedP,1);


    return particleList_[res[0].idx].iter_;
}

// =============================================================================

template<class CloudType, class mmcCloudType>
void Foam::nearestNeighbor<CloudType, mmcCloudType>::coupleSources
(
    const parcelType& p,    // Particle
    const typename parcelType::trackingData& td,
    mmcCloudType& mmcCloud,
    const scalar& mDot,
    const scalar& hDot,
    const label jFuel,
    const scalar& dt
) const
{
    auto scP = locateClosestStochasticParticle(p,td, mmcCloud);

    const labelHashTable& XiCIndexes = mmcCloud.coupling().XiC().cVarInXiC();


    if (mDot > 0.0)
    {
        // Mass and heat transfer

        // Get some properties
        const scalar m_sp = scP->m();       // mass of the stochastic particle
        const scalar Pi_sp = mDot/m_sp;     // Normalized evaporation rate

        // Mass source term
        scP->dpMsource() = mDot*dt;

        // Mixture fraction source term
        scP->dpXiCsource()[XiCIndexes["z"]] = ((1.0-scP->XiC()[XiCIndexes["z"]])*Pi_sp)*dt;

        // Species mass fraction source term
        forAll(mmcCloud.composition().componentNames(),ns)
        {
            if (ns == jFuel)
                scP->dpYsource()[ns] = ((1.0-scP->Y()[ns])*Pi_sp)*dt;
            else
                scP->dpYsource()[ns] = ((0.0-scP->Y()[ns])*Pi_sp)*dt;
        }

        // Absolute enthalpy source term
        const scalar hDqDPi = hDot/m_sp;

        scP->dp_hAsource() = (hDqDPi-scP->hA()*Pi_sp)*dt;
    }
    else
    {
        // Pure heat transfer

        // Mass of the stochastic particle
        scalar m_sp = scP->m();

        // Mass source term
        scP->dpMsource() = 0.0;

        // Mixture fraction source term
        scP->dpXiCsource()[XiCIndexes["z"]] = 0.0;

        // Species mass fraction source term
        forAll(mmcCloud.composition().componentNames(),ns)
        {
            scP->dpYsource()[ns] = 0.0;
        }

        // Absolute enthalpy
        scP->dp_hAsource() = hDot/m_sp*dt;
    }

    // Update stochastic particle properties
    scP->updateDispersedSources();
    scP->wt() = scP->m()/mmcCloud.deltaM();
    scP->T()  = mmcCloud.composition().particleMixture
    (
        scP->Y()
    ).THa(scP->hA(),scP->pc(),scP->T());

};


template<class CloudType, class mmcCloudType>
void Foam::nearestNeighbor<CloudType, mmcCloudType>::getGasPhaseConditions
(
    const parcelType& p,
    const typename parcelType::trackingData& td,
    mmcCloudType& mmcCloud,
    scalar& TInf,
    scalar& pInf,
    List<scalar>& YInf
) const
{
    auto scP = locateClosestStochasticParticle(p,td, mmcCloud);
 
    TInf = scP->T();
    pInf = scP->pc();
    YInf = scP->Y();
};

// ************************************************************************* //
