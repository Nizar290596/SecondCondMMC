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

#ifndef makePremixedMixingPopeParticleMixingModels_H
#define makePremixedMixingPopeParticleMixingModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NoMixing.H"
#include "MMCcurl.H"
#include "premixedMMCcurl.H"
#include "Curl.H"
#include "ransMMCcurl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePremixedMixingPopeParticleMixingModels(CloudType)               \
                                                                            \
    makeCloudMixingModel(CloudType);                                        \
    makeCloudMixingModelType(NoMixing, CloudType);                          \
    makeCloudMixingModelType(premixedMMCcurl, CloudType);                   \
    makeCloudMixingModelType(Curl, CloudType);                              \
    makeCloudMixingModelType(ransMMCcurl, CloudType);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
