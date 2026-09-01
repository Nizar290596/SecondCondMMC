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

#include "phiMMCcurl.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::phiMMCcurl<CloudType>::phiMMCcurl
(
    const dictionary& dict,
    CloudType& owner,
    const mmcVarSet& Xi
)
:
    MMCcurl<CloudType>(dict, owner, Xi)
{
    Info<< token::TAB << "Mixes:   phi only (composition is left to the "
        << "stage that owns it)" << endl;
}


template<class CloudType>
Foam::phiMMCcurl<CloudType>::phiMMCcurl(const phiMMCcurl<CloudType>& cm)
:
    MMCcurl<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::phiMMCcurl<CloudType>::mixPairProperties
(
    particleType& p,
    const eulerianFieldData& /* pEulerianFields */,
    particleType& q,
    const eulerianFieldData& /* qEulerianFields */,
    const scalar mixExtent,
    const scalar /* tauMix */,
    const scalar /* deltaT */
)
{
    const scalar wtSum = p.wt() + q.wt();

    if (wtSum < VSMALL)
        return;

    // Weighted pair mean relaxed by mixExtent - the same modified Curl step
    // MMCcurl applies to the composition, restricted to phi.
    const scalar phiAv = (p.wt()*p.phi() + q.wt()*q.phi())/wtSum;

    p.phi() += mixExtent*(phiAv - p.phi());
    q.phi() += mixExtent*(phiAv - q.phi());
}


// ************************************************************************* //
