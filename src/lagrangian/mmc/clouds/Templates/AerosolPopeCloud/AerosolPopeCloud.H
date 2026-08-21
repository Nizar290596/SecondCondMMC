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
    Foam::AerosolPopeCloud

Description
    Templated base class for aerosol method cloud

    - adds aerosol method
    - adds homogeneous nucleations and condensation sub-model

SourceFiles
  AerosolPopeCloudI.H
  AerosolPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef AerosolPopeCloud_H
#define AerosolPopeCloud_H

#include "ItoPopeCloud.H"
#include "aerosolPopeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declaration of classes

template<class CloudType>
class SynthesisModel;

/*---------------------------------------------------------------------------*\
Class AerosolPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AerosolPopeCloud
:
    public CloudType,
    public aerosolPopeCloud
{
    public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef AerosolPopeCloud<CloudType> aerosolPopeCloudType;
        
    private:

    // Private data

        //- Cloud copy pointer
        autoPtr<AerosolPopeCloud<CloudType> > cloudCopyPtr_;

    // Private member functions

        //- Disallow default bitwise copy construct
        AerosolPopeCloud(const AerosolPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const AerosolPopeCloud&);
        

    protected:

    // Protected data
        
        //- synthesis model
        autoPtr<SynthesisModel<AerosolPopeCloud<CloudType> > > synthesisModel_;

        // Protected Member Functions

    // Initialisation
      
        //- Set cloud sub-models
        void setModels();
        
    // Cloud evolution functions

    //- Reset state of cloud
        void cloudReset(AerosolPopeCloud<CloudType>& c);


    public:

    // Constructors

    //- Construct with explicit fields
    AerosolPopeCloud
    (
        const word& cloudName,
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& DEff,
        const volScalarField& rho,
        const volVectorField& gradRho,
        const mmcVarSet& Xi,
        const Switch initAtCnstr = true,
        bool readFields = true
    );

    //- Construct and lookup fields in object registry
    AerosolPopeCloud
    (
        const word& cloudName,
        const fvMesh& mesh,
        const mmcVarSet& Xi,
        const Switch initAtCnstr = true,
        bool readFields = true
    );


    //- Destructor
    virtual ~AerosolPopeCloud();

    // Sub-models

    //- Return const access to aerosol synthesis model
    inline const SynthesisModel<AerosolPopeCloud<CloudType> >& synthesis() const;

    //- Return access to aerosol synthesis model
    inline SynthesisModel<AerosolPopeCloud<CloudType> >& synthesis();
        
    //- Set particle thermo properties
        virtual void setParticleProperties
        (
            particleType& particle,
            const scalar& mass,
            const scalar& wt,
            const scalar& patchI,
            const scalar& patchFace,
            const bool& iniRls
        );

    //- Set Eulerian statistics capturing
    virtual void setEulerianStatistics();

    //- Update Eulerian statistics
    virtual void updateEulerianStatistics();
    
    virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AerosolPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AerosolPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
