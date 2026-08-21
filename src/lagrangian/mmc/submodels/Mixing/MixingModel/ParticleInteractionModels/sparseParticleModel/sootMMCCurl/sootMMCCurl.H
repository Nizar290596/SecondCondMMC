/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sootMMCCurl

Description
    MMC mixing model for cases with soot modeling included. Based on Curls
    model.
    
    Currently only the aISO model is supported C&K model cannot be selected

    Detail can be found in:
    [1] Son Vo, O.T. Stein, A. Kronenburg, M.J. Cleary,
        Assessment of mixing time scales for a sparse particle method,
        Combustion and Flame,
        Volume 179,
        2017,
        Pages 280-299,
        ISSN 0010-2180,
        https://doi.org/10.1016/j.combustflame.2017.02.017. 
    

\*---------------------------------------------------------------------------*/

#ifndef sootMMCCurl_H
#define sootMMCCurl_H

#include "mixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sootMMCCurl Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class sootMMCCurl
:
    public mixParticleModel<CloudType>
{
    // Typedef for particle type
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;
    
    // Private data
    
        //- Time scale model coefficients
        const scalar CL_;
        const scalar CE_;
        const scalar beta_;

        //- aISO control
        const Switch aISO_;

        //- true: harmonic mean of two particle mixing time scale
        //- false: min of two particle mixing time scale
        const Switch meanTimeScale_;
    
    public:

        //- Runtime type information
        TypeName("sootMMCCurl");

    // Constructors

        //- Construct from dictionary
        sootMMCCurl(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi);

        //- Construct copy
        sootMMCCurl(const sootMMCCurl<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new sootMMCCurl<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~sootMMCCurl() = default;
        

    // Member Functions
    
        //- Mixing particle pairs
        //  See also [1] for description of the time scale
        void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        );
        
        virtual const scalarField XiR0(label patch, label patchFace);
        virtual const scalarField XiR0(label celli);   

    // Info Statement

        //- Print an info statement to the screen
        void printInfo();   
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sootMMCCurl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
