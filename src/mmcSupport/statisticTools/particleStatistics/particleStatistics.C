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

\*---------------------------------------------------------------------------*/
#include "particleStatistics.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleStatistics::particleStatistics
(
    const dictionary& particleStatisticsDictionary,
    const fvMesh& mesh
)
    :
    mesh_(mesh),
    
    dict_(particleStatisticsDictionary), 

    particleSamplingEnabled_(dict_.lookupOrDefault<Switch>("enabled",false)),
    
    particleSamplingInterval_(dict_.lookupOrDefault<label>("sampleInterval", 1)),

    fields_(dict_.lookupOrDefault<List<word>>("fields",wordList())),

    condVarName_(dict_.lookupOrDefault<word>("condVariable","")),
    
    lowerBound_(dict_.lookupOrDefault<scalar>("lowerBound",0.0)),
    
    upperBound_(dict_.lookupOrDefault<scalar>("upperBound",1.0)),
    
    particleSampleProbesDict_
    (
        particleSamplingEnabled_ ?
        dict_.subDict("particleSampleProbes")
      : dict_ 
    ),

    collated_(dict_.lookupOrDefault<bool>("collated",true))
{
    if (particleSamplingEnabled_)
    {
        Info << "Particle Sampling enabled"
             << dict_;
    }
    setParticleSampling();
}

  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleStatistics::~particleStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particleStatistics::setParticleSampling()
{
    if (!particleSamplingEnabled_)
        return;
    
    probesName_.setSize(particleSampleProbesDict_.size());
    probes_.setSize(particleSampleProbesDict_.size());
    
    forAll(particleSampleProbesDict_.toc(),i) 
    {
        probesName_[i] = particleSampleProbesDict_.toc()[i];

        probes_[i].reset
        (
            particleProbes::New
            (
                particleSampleProbesDict_.subDict(probesName_[i])
            )
        );
    }
    
    // Set the base file path
    baseFilePath_ = mesh().time().globalPath()
        + "/mmcStatistics/particleSampling/";

    mkDir(baseFilePath_);


    // If collated IO is used the stream position must be set for a 
    // restart
    if (collated_)
    {
        // Find all files in the base file path:
        auto fileList = Foam::readDir(baseFilePath_);
        for (const fileName& fName : fileList)
        {
            if (fName.ext() == "bin")
            {
                parallelIO::setStreamPosBasedOnTime
                (
                    mesh_.time().time().value(),
                    0,
                    baseFilePath_+fName
                );
            }
        }
    }
}

// ************************************************************************* //
