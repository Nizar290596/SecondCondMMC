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

#include "mmcVarSet.H"

namespace Foam
{
    defineTypeNameAndDebug(mmcVarSet, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::mmcVarSet::mmcVarSet()
:
    varSet(),
    
    rVarIndexinXi_(),
    
    cVarIndexinXi_(),
    
    rVarIndexinXiR_(),
    
    cVarIndexinXiC_()
{}; 


Foam::mmcVarSet::mmcVarSet
(
    const dictionary& VarDict,
    const fvMesh& mesh
)
:
    varSet
    (
        VarDict,
        mesh
    ),
    
    rVarIndexinXi_(),
    
    cVarIndexinXi_(),
    
    rVarIndexinXiR_(),
    
    cVarIndexinXiC_()
{
    label cNum = 0;
    label rNum = 0;
    
    //- Initialize reference names and coupling names lists  
    forAllConstIter(varTable,this->varNames(),iter)
    {
        // Coupling variables names list
        if (isA<couplingVar>(this->Vars(*iter)))
        {
            cVarIndexinXi_.insert(this->Vars(*iter).cName(),this->varNames()[*iter]);
            cVarIndexinXiC_.insert(this->Vars(*iter).cName(),cNum);
            cNum++;
        }
        
        // Reference variables names list
        if (isA<referenceVar>(this->Vars(*iter)))
        {
            rVarIndexinXi_.insert(this->Vars(*iter).rName(),this->varNames()[*iter]);
            rVarIndexinXiR_.insert(this->Vars(*iter).rName(),rNum);
            rNum++;
        }
    }
             
    cVarIndexinXi_.resize(cVarIndexinXi_.size());
    cVarIndexinXiC_.resize(cVarIndexinXiC_.size());    
        
    rVarIndexinXi_.resize(rVarIndexinXi_.size());
    rVarIndexinXiR_.resize(rVarIndexinXiR_.size());
    
//    Info << "rVarIndexinXiR is: " << rVarIndexinXiR_ << nl
//         << "rVarIndexinXi is: " << rVarIndexinXi_   << nl
//         << "cVarIndexinXiC is: " << cVarIndexinXiC_ << nl
//         << "cVarIndexinXi is: " << cVarIndexinXi_ << endl;
}


Foam::mmcVarSet::mmcVarSet(const mmcVarSet& mmcVs)
:
    varSet(mmcVs),
    rVarIndexinXi_(mmcVs.rVarIndexinXi_),    
    cVarIndexinXi_(mmcVs.cVarIndexinXi_),
    rVarIndexinXiR_(mmcVs.rVarIndexinXiR_),    
    cVarIndexinXiC_(mmcVs.cVarIndexinXiC_)
{}; 

// ************************************************************************* //
