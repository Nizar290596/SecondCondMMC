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

#include "SinglePhaseMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SinglePhaseMixture<CloudType>::SinglePhaseMixture
(
    const dictionary& dict,
    CloudType& owner
)
:
    CompositionModel<CloudType>(owner),
    printMolFracEnabled_(dict.lookup("printMoleFractions"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SinglePhaseMixture<CloudType>::~SinglePhaseMixture()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
const Foam::wordList& Foam::SinglePhaseMixture<CloudType>::componentNames() const
{
    return this->slgThermo().carrier().species();
}


template<class CloudType>
const Foam::scalarField Foam::SinglePhaseMixture<CloudType>::YMixture0(label patch, label patchFace) 
{

    label nSpecies=this->slgThermo().carrier().Y().size();
    const PtrList<volScalarField>& fvY=this->slgThermo().carrier().Y();

    scalarField Y(nSpecies,0.0);
    for(label i=0; i<nSpecies; i++)
    {
        Y[i]=fvY[i].boundaryField()[patch][patchFace];
    }
    
    return Y;
}


template<class CloudType>
const Foam::scalarField Foam::SinglePhaseMixture<CloudType>::YMixture0(label celli) 
{

    label nSpecies=this->slgThermo().carrier().Y().size();
    const PtrList<volScalarField>& fvY=this->slgThermo().carrier().Y();

    scalarField Y(nSpecies,0.0);
    for(label i=0; i<nSpecies; i++)
    {
        Y[i]=fvY[i][celli];
    }

    return Y;
}


template<class CloudType>
const Foam::Switch Foam::SinglePhaseMixture<CloudType>::printMoleFractionsEnabled() const
{
    return printMolFracEnabled_;
}
// ************************************************************************* //

