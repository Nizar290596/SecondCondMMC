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
    Foam::PremixedMixingPopeCloud

Description
    Templated base class for mixing cloud

SourceFiles
    PremixedMixingPopeCloudI.H
    PremixedMixingPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef PremixedMixingPopeCloud_H
#define PremixedMixingPopeCloud_H

#include "ItoPopeCloud.H"
#include "ThermoPopeCloud.H"
#include "mixingPopeCloud.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class CloudMixingModel;

/*---------------------------------------------------------------------------*\
                    Class PremixedMixingPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PremixedMixingPopeCloud
:
    public CloudType,
    public mixingPopeCloud
{
    public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of particle the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef PremixedMixingPopeCloud<CloudType> premixedMixingPopeCloudType;
        
        //- Required so that the CloudMixingModel finds this type
        typedef PremixedMixingPopeCloud<CloudType> mixingPopeCloudType;

    private:

    // Private data

        //- Cloud copy pointer
        autoPtr<PremixedMixingPopeCloud<CloudType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        PremixedMixingPopeCloud(const PremixedMixingPopeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const PremixedMixingPopeCloud&);

    protected:

    // Protected data

        //- Mixing model
        autoPtr<CloudMixingModel<PremixedMixingPopeCloud<CloudType> > >  mixingModel_;

    // Protected Member Functions

        //- Set cloud sub-models
        void setModels(const mmcVarSet& Xi);

        //- Reset state of cloud
        void cloudReset(PremixedMixingPopeCloud<CloudType>& c);
        
     // Dictionaries for settings 
                    
        //- Sub-dictionary for the reference variable relaxation
        const dictionary XiRrelaxationDict_;
        
        // true: reference variable relaxation is on
        // false: reference variable relaxation is off
        const Switch XiRrelaxEnabled_; 
        
        //- Relaxation type (constantValue/dynamicRelax)
        const word XiRrelaxType_; 
        
        //- Constant relaxation factor value (when XiRrelaxFactorType_ is "constantValue")
        const scalar XiRrelaxConstValue_; 
        
      // Parameters for the dynamic relaxation
      
	//- Lower limit for the particle progress variable averaging
        const scalar XiRDynRelaxLower_;
        
        //- Upper limit for the particle progress variable averaging
        const scalar XiRDynRelaxUpper_;
        
        //- Ratio between the new and old values of the relaxation factor
        const scalar XiRDynRelaxRate_;
        
        //- Minimum number of time steps that ensures sufficient time for the flame 
        //  to adjust to new conditions prior to a further change of the relaxation factor
        const label  XiRDynStepNum_;
        
        //- Minimal value of the relaxation factor, set to avoid strong over-relaxation
        const scalar XiRDynRelaxMin_;
        
        //- Critical value of the averaged particle progress variable; 
        //  when reached the relaxation factor decreases (p.XiRrlxFactor()*XiRDynRelaxRate)
        const scalar XiRDyncPcrit_;
        
        //- Critical maximal value of the averaged particle progress variable; 
        //  when reached switch to minimal relaxation factor immediately (p.XiRrlxFactor() = XiRDynRelaxMin)
        const scalar XiRDyncPcritMax_;
        
        //- Minimal value of the averaged particle progress variable; 
        //  when reached increase relaxation factor back to 1 (p.XiRrlxFactor() = 1)
        const scalar XiRDyncPcritMin_;
        
        //- File stream to write out the relaxation parameters
        autoPtr<OFstream> dynamicRelaxOut_;

  public:

    // Constructors

    PremixedMixingPopeCloud
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

    //- Destructor
    virtual ~PremixedMixingPopeCloud();


    // Member Functions

    // Access

        //- Return a reference to the cloud copy
        inline const PremixedMixingPopeCloud& cloudCopy() const;


    // Sub-models

        //- Return const access to mixing model
        inline const CloudMixingModel<PremixedMixingPopeCloud<CloudType> >& mixing() const;

        //- Return non-const access to mixing model
        inline CloudMixingModel<PremixedMixingPopeCloud<CloudType> >& mixing();

    // I-O
    
        //- Add new particle

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
        
        //- Solve the cloud - calls all evolution functions
        template<class TrackCloudType> 
        void solve
        (
            TrackCloudType& cloud,
            typename particleType::trackingData& td
        );
        
        //- Calculate the reference variable relaxation properties
        virtual void calcXiRrelax();
};


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PremixedMixingPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PremixedMixingPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
