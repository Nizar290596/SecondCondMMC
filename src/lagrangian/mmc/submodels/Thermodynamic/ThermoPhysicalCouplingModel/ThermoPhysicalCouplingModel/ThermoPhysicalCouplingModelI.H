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
inline const Foam::boolList& Foam::ThermoPhysicalCouplingModel<CloudType>::solveEqvSpecie() const
{
    return solveEqvSpecie_;
}


template<class CloudType>
inline const Foam::mmcVarSet&
Foam::ThermoPhysicalCouplingModel<CloudType>::XiC() const
{
    return XiC_;
}



template<class CloudType>
inline const Foam::Switch& Foam::ThermoPhysicalCouplingModel<CloudType>::tauRelaxBlending() const
{
    return tauRelaxBlending_;
}


template<class CloudType>
inline const Foam::dimensionedScalar& Foam::ThermoPhysicalCouplingModel<CloudType>::tauRelaxStart() const
{
    return tauRelaxStart_;
}


template<class CloudType>
inline const Foam::dimensionedScalar& Foam::ThermoPhysicalCouplingModel<CloudType>::tauRelaxTarget() const
{
    return tauRelaxTarget_;
}


template<class CloudType>
inline const Foam::scalar& Foam::ThermoPhysicalCouplingModel<CloudType>::tauRelaxDelta() const
{
    return tauRelaxDelta_;
}


template<class CloudType>
inline const Foam::word& Foam::ThermoPhysicalCouplingModel<CloudType>::tauUnits() const
{
    return tauUnits_;
}


// New functions to aid setParticle in Clouds
template<class CloudType>
inline const Foam::scalarField Foam::ThermoPhysicalCouplingModel<CloudType>::XiC0(label patch, label patchFace)
{
    scalarField XXi(XiC_.numXiC(),0.0);

    const labelHashTable& XiIndexes = XiC_.cVarInXi();
    const labelHashTable& XiCIndexes = XiC_.cVarInXiC();

    forAllConstIter(wordList,XiCNames_,iter)
    {
        XXi[XiCIndexes[*iter]] = XiC_.Vars(XiIndexes[*iter]).field().boundaryField()[patch][patchFace];
    }
    return XXi;
}


template<class CloudType>
inline const Foam::scalarField Foam::ThermoPhysicalCouplingModel<CloudType>::XiC0(label celli)
{
    scalarField XXi(XiC_.numXiC(),0.0);

    const labelHashTable& XiIndexes = XiC_.cVarInXi();
    const labelHashTable& XiCIndexes = XiC_.cVarInXiC();

    forAllConstIter(wordList,XiCNames_,iter)
    {
        XXi[XiCIndexes[*iter]] = XiC_.Vars(XiIndexes[*iter]).field()[celli];
    }
    return XXi;
}


template<class CloudType>
inline const Foam::wordList&
Foam::ThermoPhysicalCouplingModel<CloudType>::XiCNames() const
{
    return XiCNames_;
}


template<class CloudType>
inline const Foam::word&
Foam::ThermoPhysicalCouplingModel<CloudType>::cVarName() const
{
    return cVarName_;
}


//- Return const access to Indicator
template<class CloudType>
inline const Foam::volScalarField& Foam::ThermoPhysicalCouplingModel<CloudType>::Indicator() const
{
    return Indicator_;
}


//- Return non-const access to Indicator
template<class CloudType>
inline Foam::volScalarField&
Foam::ThermoPhysicalCouplingModel<CloudType>::Indicator()
{
    return Indicator_;
}

template<class CloudType>
inline bool
Foam::ThermoPhysicalCouplingModel<CloudType>::enabled() const
{
    return true;
}
