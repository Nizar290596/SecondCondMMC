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

Class
    Foam::NoSoot

Description
    Dummy class for 'none' option - will raise an error if any functions are
    called that require return values.

SourceFiles
    NoSoot.C

\*---------------------------------------------------------------------------*/

#ifndef NoSoot_H
#define NoSoot_H

#include "SootModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class NoSoot Declaration
\*---------------------------------------------------------------------------*/
template<class CloudType>
class NoSoot
:
    public SootModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    //- Constructors

        //- Construct from dictionary
        NoSoot(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoSoot(const NoSoot<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SootModel<CloudType> > clone() const
        {
            return autoPtr<SootModel<CloudType> >
            (
                new NoSoot<CloudType>(*this)
            );
        }

    //- Destructor
    virtual ~NoSoot();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoSoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
