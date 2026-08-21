/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::sprayThermoParcel<ParcelType>::trackingData::trackingData
(
    TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),

    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    fInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.f()
        )
    ),
    PInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.p()
        )
    ),
    sqrGradInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.magSqrGradf()
        )
    )
{}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::sprayThermoParcel<ParcelType>::trackingData::TInterp() const
{
    return TInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::sprayThermoParcel<ParcelType>::trackingData::fInterp() const
{
    return fInterp_();
}

template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::sprayThermoParcel<ParcelType>::trackingData::PInterp() const
{
    return PInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::sprayThermoParcel<ParcelType>::trackingData::sqrGradInterp() const
{
    return sqrGradInterp_();
}


// ************************************************************************* //
