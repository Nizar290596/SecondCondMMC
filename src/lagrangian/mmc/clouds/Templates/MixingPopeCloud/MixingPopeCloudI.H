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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::CloudMixingModel<Foam::MixingPopeCloud<CloudType> >&
Foam::MixingPopeCloud<CloudType>::mixing() const
{
    return mixingModel_();
}


template<class CloudType>
inline Foam::CloudMixingModel<Foam::MixingPopeCloud<CloudType> >&
Foam::MixingPopeCloud<CloudType>::mixing()
{
    return mixingModel_();
}

template<class CloudType>
inline bool Foam::MixingPopeCloud<CloudType>::secondCondMixingEnabled() const
{
    return secondCondMixingModel_.valid();
}
 
 
template<class CloudType>
inline Foam::CloudMixingModel<Foam::MixingPopeCloud<CloudType> >&
Foam::MixingPopeCloud<CloudType>::secondCondMixing()
{
    return secondCondMixingModel_();
}
 
 
template<class CloudType>
inline const Foam::CloudMixingModel<Foam::MixingPopeCloud<CloudType> >&
Foam::MixingPopeCloud<CloudType>::secondCondMixing() const
{
    return secondCondMixingModel_();
}

// Second conditioning accessors
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondR() const
{
    return secondCondR_;
}
 
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondBeta() const
{
    return secondCondBeta_;
}
 
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondTauOU() const
{
    return secondCondTauOU_;
}
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondTu() const
{
    return secondCondTu_;
}
 
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondTb() const
{
    return secondCondTb_;
}
 
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondAPhi() const
{
    return secondCondAPhi_;
}
 
 
template<class CloudType>
inline Foam::scalar Foam::MixingPopeCloud<CloudType>::secondCondZPhi() const
{
    return secondCondZPhi_;
}
 
// ************************************************************************* //

