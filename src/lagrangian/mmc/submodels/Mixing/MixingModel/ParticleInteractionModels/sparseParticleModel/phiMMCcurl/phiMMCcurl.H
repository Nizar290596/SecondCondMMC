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

Class
    Foam::phiMMCcurl

Description
    First-conditioning mixing model that mixes ONLY the reaction progress
    variable phi.

    Identical to MMCcurl in every other respect - same particle list, same
    reference space, same pairing, same aISO time scale - it differs only in
    which particle properties the pair mixing acts on.

    This is the companion of secondCondMMCcurl in a two-level conditioning
    setup: every particle carries phi and has it mixed here, while the
    composition (Y, T, hA) of the second-conditioning subset is mixed by
    secondCondMMCcurl in its own reference space. Selecting it explicitly, as
    a model in its own right, is what keeps the two stages independent -
    neither model tests for the presence of the other, and each is configured
    entirely by its own coefficients sub-dictionary.

    Usage in cloudProperties:
    \verbatim
    subModels
    {
        mixingModel             phiMMCcurl;
        secondCondMixingModel   secondCondMMCcurl;

        phiMMCcurlCoeffs
        {
            r_i             1.0E-3;

            Xim_i                       // one "<axis>_m" per pairing axis
            {
                sPx_m       0.05;
                sPy_m       0.05;
                sPz_m       0.05;
            }

            pairingMethod   global;

            CL              0.5;
            CE              0.1;
            aISO            true;
            meanTimeScale   true;
        }
    }
    \endverbatim

    Note that mixing only phi here means Y, hA and the coupling variables are
    not mixed by the first conditioning at all. That is the intended division
    of labour in a two-level setup, but it does mean particles outside the
    second-conditioning subset keep their injection composition. Use plain
    MMCcurl if you want the usual single-level behaviour.

\*---------------------------------------------------------------------------*/

#ifndef phiMMCcurl_H
#define phiMMCcurl_H

#include "MMCcurl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class phiMMCcurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class phiMMCcurl
:
    public MMCcurl<CloudType>
{
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData
        eulerianFieldData;


protected:

        //- Mix the progress variable only.
        //  phi is a first-level quantity carried by every particle; the
        //  composition is left to whichever stage owns it.
        virtual void mixPairProperties
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            const scalar mixExtent,
            const scalar tauMix,
            const scalar deltaT
        ) override;


public:

    //- Runtime type information
    TypeName("phiMMCcurl");


    // Constructors

        //- Construct from dictionary
        phiMMCcurl
        (
            const dictionary& dict,
            CloudType& owner,
            const mmcVarSet& Xi
        );

        //- Construct copy
        phiMMCcurl(const phiMMCcurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType>> clone() const
        {
            return autoPtr<CloudMixingModel<CloudType>>
            (
                new phiMMCcurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~phiMMCcurl() = default;
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "phiMMCcurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
