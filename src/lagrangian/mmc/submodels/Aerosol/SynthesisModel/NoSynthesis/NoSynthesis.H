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
    Foam::NoSynthesis

Description
    Dummy class for 'none' option - will raise an error if any functions are
    called that require return values.

SourceFiles
    NoSynthesis.C

\*---------------------------------------------------------------------------*/

#ifndef NoSynthesis_H
#define NoSynthesis_H

#include "SynthesisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoSynthesis Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoSynthesis
:
    public SynthesisModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from dictionary
        NoSynthesis(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoSynthesis(const NoSynthesis<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SynthesisModel<CloudType> > clone() const
        {
            return autoPtr<SynthesisModel<CloudType> >
            (
                new NoSynthesis<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoSynthesis();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoSynthesis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
