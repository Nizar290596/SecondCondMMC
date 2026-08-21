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
    Foam::particleTrack

Description
    Gather statistics of the particles and write them into a log file 
    



SourceFiles
    particleTrack.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>
    Sergio Gutiérrez Sánchez <sergio.gutierrez-sanchez@itv.uni-stuttgart.de>

\*---------------------------------------------------------------------------*/

#ifndef particleTrack_H
#define particleTrack_H

#include <random>
#include <set>
#include "IOdictionary.H"
#include "fvMesh.H"
#include "parallelFileIO.H"
#include "particleProbes.H"
#include "particleStatistics.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                        Class particleTrack Declaration
\*---------------------------------------------------------------------------*/
class particleTrack
{
    private:
     
    // Private Member Variables

        const fvMesh& mesh_;
        
        //- Particle tracking dictionary 
        const dictionary dict_;
        
        const Switch particleTrackEnabled_;  
        const label particleTrackInterval_;

        // Number of particle to track per core
        const label particlesTrackedPerCore_;

        //- List of variables to sample
        //  if none is given, all particle data will be sampled
        const wordList fields_;

        //- Write collated 
        bool collated_;

        //- Check if a cloud is already set to track
        bool cloudSet_{false};
        
        //- Particles identification
        List<List<label>> particlesOrigID_;
        List<List<label>> particlesOrigProcID_;

        //- Path to the directory containing particle information
        fileName baseFilePath_;

    // Private Member Functions

        
        
        const fvMesh& mesh() {return mesh_;}

    public:

        // Constructors

        particleTrack
        (
            const dictionary& particleTrackDictionary,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        particleTrack(const particleTrack&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const particleTrack&) = delete;
     
        // Destructor
        virtual ~particleTrack();


        //- Set particle sampling properties
        template<class CloudType>
        void addCloudsToTrack(const CloudType& cloud);

        //- Set particle track properties
        template<class CloudType>
        void particleTracking(const CloudType& cloud);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "particleTrackTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
