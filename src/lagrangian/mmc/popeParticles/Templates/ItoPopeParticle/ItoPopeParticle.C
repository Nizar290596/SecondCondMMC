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

Class
    ItoPopeParticle

Description

\*----------------------------------------------------------------------------*/
#include "ItoPopeParticle.H"
#include "IOstreams.H"
#include "meshTools.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct as a copy
template<class ParticleType>
Foam::ItoPopeParticle<ParticleType>::ItoPopeParticle 
(
    const ItoPopeParticle<ParticleType>& p
)
:
    ParticleType(p),
    
    sCell_(p.sCell_),
    
    mass_(p.mass_),
    
    dpMsource_(p.dpMsource_),
    
    wt_(p.wt_),
    
    A_(p.A_),
    
    rw_(p.rw_),
    
    toPos_(p.toPos_)
{
    initStatisticalSampling();
}


//- Construct as a copy
template<class ParticleType>
Foam::ItoPopeParticle<ParticleType>::ItoPopeParticle
(
    const ItoPopeParticle<ParticleType>& p,
    const polyMesh& mesh
)
:
    ParticleType(p,mesh),
    
    sCell_(p.sCell_),
    
    mass_(p.mass_),
    
    dpMsource_(p.dpMsource_),
    
    wt_(p.wt_),
    
    A_(p.A_),
    
    rw_(p.rw_),
    
    toPos_(p.toPos_)
{
    initStatisticalSampling();
}
    

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    tetIndices tetIs = this->currentTetIndices();
    
    //- Super Cell
    sCell() = cloud.pManager().getSuperCellID(cellI);
    vector Upos = td.UInterp().interpolate(this->coordinates(),tetIs);

    //- Density
    scalar rho = td.rhoInterp().interpolate(this->coordinates(),tetIs);
    
    //- Grad(rho * D)
    vector gradRho = td.gradRhoInterp().interpolate(this->coordinates(),tetIs);

    //- Pope particle drift coefficient
    A() = Upos + (gradRho/rho);
    
    //- Pope particle randomwalk
    scalar rc = max(td.DInterp().interpolate(this->coordinates(),tetIs),0);

    rw() = vector
           (
                cloud.rndGen().Normal(0,1)*sqrt(dt)*sqrt(2*rc),
                cloud.rndGen().Normal(0,1)*sqrt(dt)*sqrt(2*rc),
                cloud.rndGen().Normal(0,1)*sqrt(dt)*sqrt(2*rc)
            );
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI
)
{
    this->wt() = this->m()/cloud.deltaM();
}


template<class ParticleType>
template<class TrackCloudType>
bool Foam::ItoPopeParticle<ParticleType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::particleType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::particleType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);


        ttd.switchProcessor = false;
        ttd.keepParticle = true;

        scalar f = 1.0 - p.stepFraction();
                        
        //- First order Lagrangian scheme, cell values updated only once per 
        //  time step
        
        if (p.stepFraction() == 0)
        {
            p.setCellValues(cloud,ttd,trackTime,p.cell());
        }
        
        //- toPos changed in case of cyclic BCs - should not change otherwise 
        //  for stepFractions /=0,1
        dxDeterministic() = A_ * f * trackTime;
        
        toPos() = p.position() +  dxDeterministic() + rw_* f;
        
        
        while (ttd.keepParticle && !ttd.switchProcessor && p.stepFraction() < 1)
        {
            // Cache the current position, cell and step-fraction
            const point start = p.position();
            const scalar sfrac = p.stepFraction();
            const label celli = p.cell();
     
            // Total displacement over the time-step
            // --> This differs from OpenFOAM standard KinematicParcel.C
            const vector s = toPos() - start;
     
            // Deviation from the mesh centre for reduced-D cases
            const vector d = p.deviationFromMeshCentre();

            scalar f = 1 - p.stepFraction();
            
            // --> ENHANCEMENT <--
            // Potentially add for the future
            //f = min(f, maxCo_);
            //f = min(f, maxCo_*l/max(SMALL*l, mag(s)));

            p.trackToFace(s - d, f);
            const scalar dt = (p.stepFraction() - sfrac)*trackTime;

            
            // Avoid problems with extremely small timesteps
            if (dt > ROOTVSMALL)
            {
                // Solve the chemistry on the particle
                p.calc(cloud, ttd, dt, celli);
            }
           
            if (p.onFace() && ttd.keepParticle)
            {
                p.hitFace(s, cloud, ttd);
            }
        }

	if (ttd.switchProcessor && this->cell() >= 0)
        {
            Foam::barycentric& coords =
                const_cast<Foam::barycentric&>(this->coordinates());

            // Find the smallest (should be near-zero) coordinate
            label zeroIdx = -1;
            scalar minVal = GREAT;
            for (label i = 0; i < 4; ++i)
            {
                if (coords[i] < minVal)
                {
                    minVal = coords[i];
                    zeroIdx = i;
                }
            }

            // Only nudge if a coord is actually near-zero (particle truly
            // on a face). Skip if all coords are >> 0 (particle is in
            // cell interior somehow, in which case the deadlock root cause
            // is something else and we shouldn't perturb).
            if (zeroIdx >= 0 && minVal < 1.0e-6)
            {
                const scalar nudge = 1.0e-3;
                coords[zeroIdx] = nudge;
                for (label i = 0; i < 4; ++i)
                {
                    if (i != zeroIdx) coords[i] -= nudge / 3.0;
                }
            }
        }

    return ttd.keepParticle;
}

template<class ParticleType>
template<class TrackCloudType>
bool Foam::ItoPopeParticle<ParticleType>::hitPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
 

    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];
    // Note: In case of a processor or wedge patch the sepcific 
    //       functions in particle::hitBoundaryFace() should be called.
    //       Otherwise, the reboundParticle function can lead to an infinite 
    //       loop
    if (isA<processorPolyPatch>(pp) || isA<wedgePolyPatch>(pp))
    {
        // Skip processor patches
        return false;
    }
    else if (pp.physicalType() == "noParticleOutflow")
    {
        toPos() = reboundParticle(pp,this->toPos(),this->position());

        td.countRebound();
        
        td.keepParticle = true;
        
        return true;
    }
    else if (pp.physicalType() == "particleInflowOutflow")
    { 
        const ParticleType& p =
            static_cast<const ParticleType&>(*this);

        //- Wall normal vector
        vector nw = pp.faceAreas()[pp.whichFace(p.face())];
    
        //- Unit wall normal
        nw /= mag(nw);

        //- Wall normal velocity
        tetIndices tetIs = this->currentTetIndices();
    
        scalar velN = td.UInterp().interpolate(this->coordinates(),tetIs) & nw;

        //- Rebound if there is inflow at the boundary
        if (velN < 0)
        {
            toPos() = reboundParticle(pp,this->toPos(),this->position());

            td.countRebound();
        }
        else
        {
            const vector toPosTmp = toPos();

            toPos() = noRandomOutflowRebound
                        (
                            pp,
                            this->toPos(),
                            this->position(),
                            dxDeterministic()
                        );

            if (toPos() == toPosTmp)
            {
                td.keepParticle = false;

                td.countLeaving();
            }
            else
                td.countNoRandomOutflow();
        }
    }
    else if (pp.physicalType() == "particleInflowFreeOutflow")
    {
        
        const ParticleType& p =
            static_cast<const ParticleType&>(*this);

        //- Wall normal vector
        vector nw = pp.faceAreas()[pp.whichFace(p.face())];
    
        //- Unit wall normal
        nw /= mag(nw);

        //- Wall normal velocity
        tetIndices tetIs = this->currentTetIndices();
    
        scalar velN = td.UInterp().interpolate(this->coordinates(),tetIs) & nw;

        //- Rebound if there is inflow at the boundary
        if (velN < 0)
        {
            toPos() = reboundParticle(pp,this->toPos(),this->position());

            td.countRebound();
        }
        else
        {
            td.keepParticle = false;

            td.countLeaving();
        }
    }
    else if (pp.physicalType() == "noRandomOutflow")
    {
        const vector toPosTmp = toPos();

        toPos() = noRandomOutflowRebound
                    (
                        pp,
                        this->toPos(),
                        this->position(),
                        dxDeterministic()
                    );

        if (toPos() == toPosTmp)
        {
            td.keepParticle = false;

            td.countLeaving();
        }
        else
            td.countNoRandomOutflow();
    }
    else if (pp.physicalType() == "particleOutflow")
    {       
        td.keepParticle = false;
    
        td.countLeaving();
    }
    else
    {
        // If it is none of these types loop over all possible patch
        // types in the Foam::particle::hitFace() function
        return false;
    }

    return true;
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::hitCyclicPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
 
    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    if (pp.physicalType() == "periodic")
    {
        
        //- Old position at cyclic boundary
        const point posOldBound = this->position();
        
        ParticleType::hitCyclicPatch(static_cast<const cyclicPolyPatch&> (pp), td);

        toPos() = periodicParticle(this->toPos(),this->position(),posOldBound);

        td.countPeriodic();
    }
    else
    {
        FatalErrorIn
        (
            "ItoPopeParticle<ParticleType>::hitCyclicPatch"
        )   << "Undefined physicalType() for cyclic patch." << nl
            << "Valid physicalTypes in blockMesh are : periodic." << nl
            << exit(FatalError); 
    }
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
 
    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    toPos() = reboundParticle(pp,toPos(),this->position());

    td.countWall();        
} 


template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::hitProcessorPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    td.switchProcessor=true;
    
    td.countChangedProzessor();
} 


template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::hitWedgePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    FatalErrorInFunction
        << "Hitting a wedge patch is not supported."
        << abort(FatalError);


    // Potential Idea
    // ==============
    // Particles that hit the wegde patch behave similar to a slip condition
    // of a wall. 
    // The wall normal component gets removed and particles can move in the
    // wall parallel directio
    // This however still leads to infinite loops for some particles, even if
    // they do not hit a wedge patch

    // typename TrackCloudType::parcelType& p =
    //     static_cast<typename TrackCloudType::parcelType&>(*this);

    // const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    // // Wall normal vector
    // vector nw = pp.faceAreas()[pp.whichFace(p.face())];
    
    // // Unit wall normal
    // nw /= mag(nw);

    // const Foam::point toPos = this->toPos();
    // const point& pos = this->position();

    // // Subtract the wall normal component of the toPos value
    // // Normal component of toPos
    // scalar toPosN = (this->toPos() - pos) & nw; //Scalar projection 

    // // Add a tiny fraction to move the particle away from the wall.
    // A_ = A_ - 1.01*(A_ & nw) * nw;
    // rw_ = rw_ - 1.01*(rw_ & nw) * nw;
    // this->toPos() = this->toPos() - 1.01*toPosN*nw;

    // td.countWedge();
}


// symmetryPlane Template added by AV
template<class ParticleType>
template<class TrackCloudType>
void Foam::ItoPopeParticle<ParticleType>::hitSymmetryPlanePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    toPos() = reboundParticle(pp,toPos(),this->position());

    td.countSymmetryPlane();
}



template<class ParticleType>
void Foam::ItoPopeParticle<ParticleType>::transformProperties (const tensor& T)
{
    ParticleType::transformProperties(T);
    
    A_ = transform(T, A_);
}


template<class ParticleType>
void Foam::ItoPopeParticle<ParticleType>::transformProperties(const vector& separation)
{
    ParticleType::transformProperties(separation);   
}


template<class ParticleType>
Foam::point Foam::ItoPopeParticle<ParticleType>::reboundParticle
(
    const polyPatch& pp,
    const Foam::point& toPos,
    const Foam::point& pos
) 
{
    const ParticleType& p =
        static_cast<const ParticleType&>(*this);

    // Wall normal vector
    vector nw = pp.faceAreas()[pp.whichFace(p.face())];
    
    // Unit wall normal
    nw /= mag(nw);

    // Normal component of toPos
    scalar toPosN = (toPos - pos) & nw; //Scalar projection 

    // Note: If phe particle is at the conjunction of a 
    // processor boundary and a rebound boundary particles can get
    // stuck because toPos is newly calculated after transferring 
    // the particle
 
    if (toPosN > 0)
    {
        A_ = A_ - 2*(A_ & nw) * nw;
        rw_ = rw_ - 2*(rw_ & nw) * nw;

        return toPos - 2*toPosN*nw;
    }
    else
    {   
        return toPos;
    }

}


template<class ParticleType>
Foam::point Foam::ItoPopeParticle<ParticleType>::noRandomOutflowRebound
(
    const polyPatch& pp,
    const Foam::point& toPos,
    const Foam::point& pos,
    const Foam::vector& dxDet
) const
{
    const ParticleType& p =
        static_cast<const ParticleType&>(*this);

    //- Wall normal vector
    vector nw = pp.faceAreas()[pp.whichFace(p.face())];
    
    //- Unit wall normal
    nw /= mag(nw);


    scalar toPosDetN = dxDet & nw;
    
    // if the deterministic walk is inside the domain, only use deterministic
    // otherwise keep the original position
    if (toPosDetN < 0 )
        return pos + dxDet;
    else
        return toPos;
}


template<class ParticleType>
Foam::point Foam::ItoPopeParticle<ParticleType>::periodicParticle
(
    const Foam::point& toPos,
    const Foam::point& pos,         //new position on opposite face
    const Foam::point& posOldBound  //position of patch hit - position before OF hitCyclicPatch
) const
{
    return toPos + (pos-posOldBound);
}


template<class ParticleType>
template<class CloudType>
void Foam::ItoPopeParticle<ParticleType>::setStaticProperties(CloudType& c)
{
   // Do nothing
}


template<class ParticleType>
void Foam::ItoPopeParticle<ParticleType>::initStatisticalSampling()
{
    nameVariableLookUpTable().addNamedVariable("mass",mass_);
    nameVariableLookUpTable().addNamedVariable("weight",wt_);
}


template<class ParticleType>
Foam::DynamicList<Foam::scalar> 
Foam::ItoPopeParticle<ParticleType>::getStatisticalData
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
Foam::ItoPopeParticle<ParticleType>::getStatisticalDataNames
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

#include "ItoPopeParticleIO.C"

// ************************************************************************* //
