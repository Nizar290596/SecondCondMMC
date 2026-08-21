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

template<class CloudType>
void Foam::particleTrack::addCloudsToTrack(const CloudType& cloud)
{
    if (!particleTrackEnabled_)
        return;
    
    if (!cloudSet_)
        cloudSet_=true;
    else
        FatalError 
            << "Only one cloud can be selected for particle tracking"
            << exit(FatalError);

    particlesOrigID_.setSize(Pstream::nProcs());
    particlesOrigProcID_.setSize(Pstream::nProcs());

    for (label i=0; i < Pstream::nProcs(); i++)
    {
        particlesOrigID_[i].setSize(particlesTrackedPerCore_);
        particlesOrigProcID_[i].setSize(particlesTrackedPerCore_);
    }

    const label N_Particles = cloud.size();
    Info << "Number of particles per core subcloud: " << N_Particles << endl;
    std::random_device rdEng; // Define a random number from hardware
    std::mt19937 eng(rdEng()); // Seed the random number generator

    std::uniform_int_distribution<> distr(0, N_Particles); // Define the range for the random numbers generated
    int k = Pstream::myProcNo();
    for (int i=0; i < particlesTrackedPerCore_; i++)
    {
        particlesOrigID_[k][i] = distr(eng);
        particlesOrigProcID_[k][i] = k;
    }
   
    Pstream::gatherList(particlesOrigID_);
    Pstream::scatterList(particlesOrigID_);

    Pstream::gatherList(particlesOrigProcID_);
    Pstream::scatterList(particlesOrigProcID_);

    Info << "List of particle indexes to do the data mining after merging: " << endl;
    Info << particlesOrigID_ << endl;

    Info << "List of particle processors indexes to do the data mining after merging: " << endl;
    Info << particlesOrigProcID_ << endl;

    // Set the base file path for particle tracking
    baseFilePath_ = mesh().time().globalPath()
        + "/mmcStatistics/particleTrack/";

    mkDir(baseFilePath_);
}

template<class CloudType>
void Foam::particleTrack::particleTracking(const CloudType& cloud)
{
    const label timeStep = mesh().time().timeIndex();

    if (!particleTrackEnabled_ || (timeStep % particleTrackInterval_))
        return;
    
    Info << "Sample particle track data..." << endl;

    DynamicList<scalarList> particlesTracked;
    for (auto& p : cloud)
    {
        label origID = p.origId();
        label origProc = p.origProc();

        forAll(particlesOrigProcID_, k)
        {
            if (origProc != k)
                continue;
            else
            {
                forAll(particlesOrigID_[k], i)
                {
                    if (origID == particlesOrigID_[k][i])
                        particlesTracked.append(p.getStatisticalData(fields_));

                }
            }
        }
    }

    List<word> header = cloud.begin()->getStatisticalDataNames(fields_);

    if (collated_)
    {
        parallelIO::writeFileBinary
        (
            particlesTracked,
            header,
            baseFilePath_ + cloud.name() + "_" + "TrackingData"
            + ".bin"
        );
    }
    else
    {
        // Write out the file 
        parallelIO::writeFile
        (
            particlesTracked,
            header,
            baseFilePath_ + cloud.name() + "_" + "TrackingData"
            + "_" + name(mesh_.time().timeIndex()) + ".dat"
        );
    }

}

// ************************************************************************* //
