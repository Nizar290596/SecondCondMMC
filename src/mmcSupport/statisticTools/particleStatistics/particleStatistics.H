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
    Foam::particleStatistics

Description
    Gather statistics of the particles and write them into a log file 
    
    All settings are read from a dictionary

    particleSampling
    {
        enabled        true;    // enable or disable particle sampling
        condVariable      z;    // Optional reference variable to limit 
                                // particle selection
        lowerBound      0.0;    // Lower limit of the reference variable
        upperBound      1.0;    // Upper limit of the reference variable

        collated        true;   // Write the data in a collated binary format
                                // default is true

        // Optionally define a set of variable fields to sample.
        fields          (rho T dx);

        // Define in this sub dictionary all probes where particles
        // are sampled. Multiple probes can be defined
        particleSampleProbes
        {
            probe1
            {
                type    particleProbeBox;   // Type of probe
                p1      (0 0 0);            // particleProbeBox specific
                p2      (1 1 1);            // particleProbeBox specific
            }

            probe2
            {
                type    particleProbeBox;   // Type of probe
                p1      (0 0 0);            // particleProbeBox specific
                p2      (1 1 1);            // particleProbeBox specific
            }

        }
    }

SourceFiles
    particleStatistics.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#ifndef particleStatistics_H
#define particleStatistics_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "parallelFileIO.H"
#include "particleProbes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                        Class particleStatistics Declaration
\*---------------------------------------------------------------------------*/
class particleStatistics
{
    private:
     
    // Private Member Variables

        const fvMesh& mesh_;
        
        //- Particle sampling dictionary 
        const dictionary dict_;
        
        const Switch particleSamplingEnabled_;  
        const label particleSamplingInterval_;

        //- List of variables to sample
        //  if none is given, all particle data will be sampled
        const wordList fields_;

        //- Name of the conditioning variable
        const word condVarName_;

        //- Lower and upper bound of the optional conditioning variable
        //  This is used to limit the particle selection for sampling
        const scalar lowerBound_;
        const scalar upperBound_;

        //- Dictionary for the particle probes
        const dictionary particleSampleProbesDict_;        

        //- Write collated 
        bool collated_;
        
        //- Name of the probes
        List<word> probesName_;
        
        //- List of particle probes 
        List<autoPtr<particleProbes>> probes_;
        
        //- Path to the directory containing particle information
        fileName baseFilePath_;

    // Private Member Functions

        //- Set particle sampling properties
        virtual void setParticleSampling();
        
        const fvMesh& mesh() {return mesh_;}
    
    public:

        // Constructors

        particleStatistics
        (
            const dictionary& particleStatisticsDictionary,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        particleStatistics(const particleStatistics&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const particleStatistics&) = delete;
     
        // Destructor
        virtual ~particleStatistics();

        //- Set particle sampling properties
        template<class CloudType>
        void particleSampling(const CloudType& cloud);
            
};

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "particleStatisticsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
