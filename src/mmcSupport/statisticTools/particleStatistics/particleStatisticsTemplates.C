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
void Foam::particleStatistics::particleSampling(const CloudType& cloud)
{
    const label timeStep = mesh().time().timeIndex();

    if (!particleSamplingEnabled_ || (timeStep % particleSamplingInterval_))
        return;
    
    Info << "Sample particle data..." << endl;

    // Loop over the probes and check if the particle is located in the probe
    forAll(probes_,i)
    {
        const particleProbes& probe = probes_[i]();
        
        DynamicList<scalarList> particlesInProbe;
        
        if (condVarName_ != "")
        {
            for (auto& p : cloud)
            {
                auto nameVarTable = p.nameVariableLookUpTable();

                scalar condVarValue = nameVarTable.get(condVarName_);

                if 
                (
                    condVarValue > lowerBound_ 
                 && condVarValue < upperBound_ 
                 && probe.containsParticle(p)
                )
                    particlesInProbe.append(p.getStatisticalData(fields_));
            }
        }
        else
        {
            for (auto& p : cloud)
            {
                particlesInProbe.append(p.getStatisticalData(fields_));
            }
        }
        

        List<word> header;

        if (cloud.size() > 0)
        {
            header = cloud.begin()->getStatisticalDataNames(fields_);
        }

        // Now it gets tricky. If there are no particles in the cloud
        // on the processor the getStatisticalDataNames cannot be called. 
        // Try to collect all the headers from all processors
        List<List<word>> allProcData(Pstream::nProcs());
        allProcData[Pstream::myProcNo()] = header;
        Pstream::gatherList(allProcData);
        Pstream::scatterList(allProcData);

        // Loop over list and make sure that it either has 0 entries in case 
        // of no particles on that processor, or n entries consistently
        for (auto& procData : allProcData)
        {
            if (procData.size() > 0)
            {
                header = procData;
                break;
            }
        }



        if (collated_)
        {
            parallelIO::writeFileBinary
            (
                particlesInProbe,
                header,
                baseFilePath_ + cloud.name() + "_" 
              + probesName_[i] + ".bin"
            );
        }
        else
        {
            // Write out the file 
            parallelIO::writeFile
            (
                particlesInProbe,
                header,
                baseFilePath_ + cloud.name() + "_" 
              + probesName_[i] + "_" + name(mesh_.time().timeIndex()) + ".dat"
            );
        }
    }
}

// ************************************************************************* //
