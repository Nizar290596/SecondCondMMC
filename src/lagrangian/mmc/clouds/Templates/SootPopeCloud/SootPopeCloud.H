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
    Foam::SootPopeCloud

Description
    Templated base class to include soot modeling

SourceFiles
    SootPopeCloudI.H
    SootPopeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef SootPopeCloud_H
#define SootPopeCloud_H

#include "ItoPopeCloud.H"
#include "sootPopeCloud.H"
#include "SootModel.H"
#include "RadiationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes

template<class CloudType>
class SootModel;

template<class CloudType>
class CloudRadiationModel;

/*---------------------------------------------------------------------------*\
                    Class ThermoPopeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SootPopeCloud
:
    public CloudType,
    public sootPopeCloud
{
    public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of particle the cloud was instantiated for
        typedef typename CloudType::particleType particleType;

        //- Convenience typedef for this cloud type
        typedef SootPopeCloud<CloudType> SootPopeCloudType;

    protected:

        //- Soot model 
        autoPtr<SootModel<SootPopeCloudType>> sootModel_;

        //- Radiation sub-model
        autoPtr<CloudRadiationModel<SootPopeCloudType > > radiationModel_;

        //- Include radiation
        Switch radiation_;

        //- Soot enabled, density and minimum molecular weight
        Switch sootingFlame_;
        scalar sootSolidDensity_;

        //- MW of the smallest soot, multiple definition
        List<scalar> sootMinMolWt_;

        //- size of sootMinMolWt_
        label numfVsoot_;

        //- min Kn number for "differential mixing". 
        //- the model is now very simple.
        scalar sootMinKn_;

        //- Threshold on fVsoot for intermittency []
        scalar IntermittencyfV_;

        //- Get the index of O2 and N2 in mechanism
        //  ->[0] index of O2
        //  ->[1] index of N2
        List<label> indexO2N2_;

    // Initialisation

        //- Set cloud sub-models
        void setModels();

        //- Reset state of cloud
        void cloudReset(SootPopeCloudType& c);

        //- Initialise the soot data
        void setSootData();

    public:

    // Constructors

        //- Construct with given fields
        SootPopeCloud
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

        //- Construct and lookup required fields in the object registry
        SootPopeCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const mmcVarSet& Xi,
            const Switch initAtCnstr = true,
            bool readFields = true
        );

        //- Destructor
        virtual ~SootPopeCloud();

        //- Disallow default bitwise copy construct
        SootPopeCloud(const SootPopeCloudType&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const SootPopeCloudType&) = delete;


    // Member Functions

    // Access

        //- Return const access to soot data
        inline const Switch& sootingFlame() const;
        inline const scalar& sootSolidDensity() const;
        inline const List<scalar>& sootMinMolWt() const;
        inline const label& numfVsoot() const;
        inline const scalar& sootMinKn() const;
        inline const scalar& IntermittencyfV() const;
     
        //- Return const access to the indexO2N2_ list
        inline const List<label>& indexO2N2() const;



    // Sub-models

        //- Return const access to soot model
        inline const SootModel<SootPopeCloudType>& soot() const;

        //- Return const access to radiation model
        inline const CloudRadiationModel<SootPopeCloudType>& radiation() const;
        
        //- Return non-const access to soot model
        inline SootModel<SootPopeCloudType>& soot();

        //- Return non-const access to radiation model
        inline CloudRadiationModel<SootPopeCloudType>& radiation();

    // Check

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
        

        //- Update the soot volume fraction
        scalarField calcfVsoot
        (
            const scalarField& Y,
            const scalar& T,
            const scalar& pc
        );

        //- Update the Intermittency
        scalarField calcIntermittency
        (
            const scalarField& fVsoot
        );

        //- calculate the XO2/XN2
        void calcXO2XN2
        (
            const scalarField& Y,
            scalar& XO2N2
        ) const;

        //- Radiation flag
        inline bool radiationFlag() const;


    // I-O

        virtual void writeFields() const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SootPopeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SootPopeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
