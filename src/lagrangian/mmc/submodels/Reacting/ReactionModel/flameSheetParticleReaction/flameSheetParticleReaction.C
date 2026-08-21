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

#include "flameSheetParticleReaction.H"

#include "vector.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::flameSheetParticleReaction<CloudType>::flameSheetParticleReaction
(
    const dictionary& dict,
    CloudType& owner
)
:
    ReactionModel<CloudType>(dict,owner,typeName),
    
    fSLU_(this->coeffDict().lookup("fSLU")),
    
    Tquench_(this->coeffDict().template get<scalar>("Tquench")),
    
    Y0S1_(3)

{
    //- Read in and store flamesheet species parameters
    forAll(Y0S1_, i)
    {
       Y0S1_[i].setSize(this->owner().composition().componentNames().size());
    }

    forAll(this->owner().composition().componentNames(), ns)
    {
        if (this->coeffDict().found(this->owner().composition().componentNames()[ns]))
        {
            List<scalar> Y0S1new(this->coeffDict().lookup(this->owner().composition().componentNames()[ns]));

            forAll(Y0S1new, i)
            {
                Y0S1_[i][ns] =  Y0S1new[i];
            }
        }
        else
        {
            for (label i=0; i < 3; i++)
            {
                Y0S1_[i][ns] = 0;
            }
        }
        
    }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::flameSheetParticleReaction<CloudType>::~flameSheetParticleReaction()
{}

// ************************************************************************* //

// It could be considered a renormalization to be able to use with 
// partially premixed flames. Z = z/zMax !!!

template<class CloudType>
void Foam::flameSheetParticleReaction<CloudType>::ignite()
{
    Info << "flameSheetParticleChemistry ignite" << endl;

    scalar fstoic = fSLU_[0];
            
    forAllIters(this->owner(), iter)
    {
        typename CloudType::particleType& p = iter();
        
        forAll(p.Y(), ns)
        {
            scalar pY0 = Y0S1_[0][ns];
            scalar pYstoic = Y0S1_[1][ns];
            scalar pY1 = Y0S1_[2][ns];
            
            if (p.XiC().first() < fstoic)// Still using first coupling variable
            {
                scalar K = (pYstoic - pY0) / fstoic;
    
                p.Y()[ns] = pY0 + K * p.XiC().first();//
            }
            else
            {
                scalar K = (pY1 - pYstoic) / (1.0 - fstoic);
        
                p.Y()[ns] = pYstoic + K * (p.XiC().first() - fstoic);//
            }
        }
    }    

}


template<class CloudType>
void Foam::flameSheetParticleReaction<CloudType>::calculate
(
    const scalar t0, 
    const scalar dt,
    const scalar ha, 
    const scalar pc, 
    const scalar pT,
    const scalar pz,
    scalarField& pY
)
{
    scalar fstoic = fSLU_[0];
    scalar fLFL = fSLU_[1];
    scalar fUFL = fSLU_[2];
            
    forAll(pY, ns)
    {
    
        if (pT > Tquench_)
        {
            scalar pY0 = Y0S1_[0][ns];
            scalar pYstoic = Y0S1_[1][ns];
            scalar pY1 = Y0S1_[2][ns];
            
            if ( (pz < fstoic) & (pz >= fLFL)  )
            {
                scalar K = (pYstoic - pY0) / fstoic;
        
                pY[ns] = pY0 + K * pz;
            }
            else if ( (pz >= fstoic) & (pz <= fUFL) )
            {
                scalar K = (pY1 - pYstoic) / (1.0 - fstoic);
            
                pY[ns] = pYstoic + K * (pz - fstoic);
            }
        }
    }
}


// ************************************************************************* //

