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

template<class ParticleType>
template<class TrackCloudType>
inline Foam::ItoPopeParticle<ParticleType>::trackingData::trackingData
(
    const TrackCloudType& cloud
)
:
    ParticleType::trackingData(cloud),    
    
    wallCollisions_(0), //Added by AV

    wedgeCollisions_(0), //Added by AV

    symmetryPlaneCollisions_(0),
    
    leavingModel_(0),
    
    injectedInModel_(0),
    
    changedProzessor_(0),
    
    noOutflowRebound_(0),

    noRandomOutflowRebound_(0),

    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U()
        )
    ),    
    
    DInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.DEff()
        )
    ),

    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho()
        )
    ), 
   
    gradRhoInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.gradRho()
        )
    )
{}


template<class ParticleType>
inline void Foam::ItoPopeParticle<ParticleType>::trackingData::reportCounters(label& particles)
{
    reduce(particles, sumOp<label>());
    
    reduce(wallCollisions_, sumOp<label>());

    reduce(wedgeCollisions_, sumOp<label>()); //Added by AV

    reduce(symmetryPlaneCollisions_, sumOp<label>()); //Added by AV
    
    reduce(leavingModel_, sumOp<label>());
    
    reduce(changedProzessor_, sumOp<label>());
    
    reduce(noOutflowRebound_, sumOp<label>());
    
    reduce(noRandomOutflowRebound_, sumOp<label>());
        
    Info << "\t"
         << particles               << " particles moved. " 
         << wallCollisions_         << " wall hits. " 
         << noOutflowRebound_       << " rebounds off noOutflow patches. "
         << noRandomOutflowRebound_ << " rebounds off noRandomOutflow patches."
         << nl << "\t"
         << leavingModel_           << " particles left the domain. " 
         << periodicModel_          << " particles hit periodic BCs. "
         << changedProzessor_       << " particles passed between processors. " 
         << endl;

    Info << "\t" << wedgeCollisions_ << " wedge hits" << endl;  //Added by AV

    Info << "\t" << symmetryPlaneCollisions_ << " symmetry Plane hits" << endl;  //Added by AV
    
}


// ************************************************************************* //
