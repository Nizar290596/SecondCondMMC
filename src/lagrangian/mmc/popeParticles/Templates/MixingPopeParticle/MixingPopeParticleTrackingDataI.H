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
inline Foam::MixingPopeParticle<ParticleType>::trackingData::trackingData
(
    TrackCloudType& cloud
)
:
    ParticleType::trackingData(cloud),
    XiRInterp_()
{

    if (cloud.mixing().numXiR())
    {           
        const mmcVarSet& setOfXi(cloud.mixing().XiR());        
        const labelHashTable& XiIndexes = setOfXi.rVarInXi();
        const labelHashTable& XiRIndexes = setOfXi.rVarInXiR();
        
        this->XiRInterp_.setSize(cloud.mixing().numXiR());

        forAllConstIter(wordList, cloud.mixing().XiRNames(), iter)
        {

            if (setOfXi.Vars(XiIndexes[*iter]).refType()=="interpolated")
            {                   
                XiRInterp_.set
                (
                    XiRIndexes[*iter],
                    interpolation<scalar>::New
                    (
                        cloud.solution().interpolationSchemes(),
                        setOfXi.Vars(XiIndexes[*iter]).field()
                    )
                );
            }
        }
    }    
}


template<class ParticleType>
inline const Foam::PtrList<Foam::interpolation<Foam::scalar> >&
Foam::MixingPopeParticle<ParticleType>::trackingData::XiRInterp() const
{
    return XiRInterp_;
}
// ************************************************************************* //
