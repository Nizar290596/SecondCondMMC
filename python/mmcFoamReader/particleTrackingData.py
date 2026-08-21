import numpy as np
import glob
import warnings
import os
from scipy import stats
import sys

class particleTrackingData:
    
    """Check the file type for ascii or binary"""
    def _checkFileType(self,fileName):
        fileEnding = fileName.split(".")[-1]
        if (fileEnding == "dat"):
            self._fileType = "ascii"
        elif (fileEnding == "bin"):
            self._fileType = "binary"

    """Find all time files for the given probe name a.k.a. fileName"""
    def _findFiles(self,dirPath,fileName):
        fileList = glob.glob(os.path.join(dirPath,fileName)+'*')
        return sorted(fileList)
    
    """Find all variables defined in the header"""
    def _setHeader(self,file):
        if (self._fileType == "ascii"):
            with open(file) as f:
                headerLine = f.readline()
                headerLine = headerLine.rstrip()
                # Split header line at comma position
                self._fields = headerLine.split('\t')
        if (self._fileType == "binary"):
            self._headerBytes = 0
            with open(file, mode="rb") as f:
                # Find beginning of binary block
                while (True):
                    char = f.read(1)
                    self._headerBytes = self._headerBytes + 1
                    if (char == b'\n'):
                        break
            with open(file, mode="rb") as f:
                buffer = (f.read(self._headerBytes)).decode("utf-8")
                buffer = buffer.strip()
                buffer.rstrip()

                self._fields = buffer.split('\t')
                
    """Find the number of particles sampled for each time index"""
    def _numberOfParticlesPerFile(self,file):
        numParticles = 0
        with open(file, "rb") as f:
            numParticles = sum(1 for _ in f) - 1

        return numParticles
    
    def __init__(self,dirPath,fileName):
        """
        Create an object of class particleTrackingData
        Load the particle tracking data into the class and provide functions
        to post-process and access the data
        
        Parameters
        ----------
        dirPath :  string
            Provide the path to the direcotry where the particle data is stored
        fileName : string
            Name of the particle statistics probe without the time index
            If the exact name is given only this file is processed
        Examples
        --------
        Load all time points of a particle data set
        >>> from mmcFoamReader import particleTrackingData
        >>> pData = particleTrackingData(
        >>>     "./mmcStatistics/particleSampling/",
        >>>     "probeName")
        Load only a specific time point of a probe
        >>> from mmcFoamReader import particleTrackingData
        >>> pData = particleTrackingData(
        >>>     "./mmcStatistics/particleSampling/",
        >>>     "probeName_100.dat")
        """
        self._fields = []
        self._data = np.array([])
        self._fileList = self._findFiles(dirPath,fileName)
        self._headerBytes = 0

        # Set to either ascii or binary
        self._fileType = "unknown"

        # Check if it is an ASCII or binary file
        # if binary it has the file ending .bin for ascii .dat
        self._checkFileType(self._fileList[0])

        if (self._fileType == "ascii"):
        
            # Get the number of particles stored in the first file
            self._maxIndex = 0
            self._minIndex = 0
            self._maxNumParticles = 0
            self._minNumParticles = sys.maxsize
            for f, file in enumerate(self._fileList):
                nParticles = self._numberOfParticlesPerFile(file)
                if (nParticles >= self._maxNumParticles):
                    self._maxNumParticles = nParticles
                    self._maxIndex = f
                if (nParticles <= self._minNumParticles):
                    self._minNumParticles = nParticles
                    self._minIndex = f
                
            # Get the number of time steps of the time series
            self._totalTimeSteps = len(self._fileList)
            
            # Read first entry to get header
            self._setHeader(self._fileList[0])
            
            # Pre allocate the space for data 
            self._data = np.empty((self._totalTimeSteps, self._maxNumParticles, len(self._fields)))

            # Read the data of the first time series step file
            currIndex = 0
            for f in self._fileList:
                print("Processing file : ",f)
                data = np.genfromtxt(f, dtype=float, delimiter='\t', skip_header=1, usecols=range(len(self._fields)), filling_values=np.nan)
                self._data[currIndex,:len(data),:] = data
                currIndex = currIndex+1
            
            print("Total time steps of the time series: ", self._totalTimeSteps)
            
            # Sort the data of each time step based on the origId
            self._sortedData = np.full((self._totalTimeSteps, self._minNumParticles, len(self._fields)), np.nan)
            ind, n = self.getfieldIndex("origId")
            origId_index = ind[n]
            
            ind, n = self.getfieldIndex("origProc")
            origProc_index = ind[n]

            referenceOrigId = self._data[self._minIndex,:self._minNumParticles,origId_index]
            referenceOrigProc = self._data[self._minIndex,:self._minNumParticles,origProc_index]
            
            for i in range(len(self._data[:,0,0])):
                filteredData = self.filterTimeStepData(data=self._data[i,:,:], refOrigId=referenceOrigId, refOrigProc=referenceOrigProc)
                self._sortedData[i,:,:] = filteredData
            
            ind, n = self.getfieldIndex("time")
            self._time_index = ind[n]
            
            
        '''
        if (self._fileType == "binary"):
            
            # Read first entry to get header
            self._setHeader(self._fileList[0])
            
            # Allocate space in the array
            numBytesInFile = os.path.getsize(self._fileList[0])
            
            self.numParticles = int((numBytesInFile - self._headerBytes)/(len(self._fields)*8))

            self._data = np.zeros((self.numParticles,len(self._fields)))

            # Read in the binary data
            rowIndex = 0
            with open(self._fileList[0], mode="rb") as f:
                # Find beginning of binary block
                f.read(self._headerBytes)

                # Read first data block
                while (True):
                    row = f.read(len(self._fields)*8)
                    if (len(row) < len(self._fields)*8):
                        break
                    self._data[rowIndex,:] = np.frombuffer(row,dtype='double')
                    rowIndex = rowIndex + 1
        '''
        
    def filterTimeStepData(self, data, refOrigId, refOrigProc):
        """
        Filters the data of the different time steps read
        Checks if a certain combination of origId and origProc is in the data
        provided
        It gets for the whole time series the particles of the most
        restrictive time step available (the one with the fewest particles)
        USAGE:
            data: Data of all the particles tracked for a single time step
            refOrigId: the origId selected to track along the time steps
            refOrigProc: the corresponding origProc to track along the time steps
        """
        
        ind, n = self.getfieldIndex("origId")
        origId_index = ind[n]
        
        ind, n = self.getfieldIndex("origProc")
        origProc_index = ind[n]

        count = 0
        filteredData = np.full((len(refOrigId), len(self._fields)), np.nan)
        for origId, procId in zip(refOrigId, refOrigProc):
            matches = ((data[:, origId_index] == origId) & (data[:, origProc_index] == procId))
            matching_indices = np.where(matches==True)[0]
            if (matching_indices.size != 0):    
                filteredData[count,:] = data[matching_indices[0],:]
            count = count+1
        
        return filteredData
    
    
    def getfieldIndex(self,name,n=-1):
        """
        Get the field data index for the given name
        USAGE:
            Returns first occurance of the fieldname
            getfieldIndex("nameOfField")
            If multiple occurances with the same name exist the 
            parameter n has to be given to specify which variable should be 
            returned.
            Note, n is zero based, e.g., getFieldIndex("z",1) returns the 
            conditioning variable not the physical z coordinate 
        """
        # get the index in the header list
        ind = []
        i = 0
        for e in self._fields:
            if (e == name):
                ind.append(i)
            i = i+1
        
        if len(ind) == 0:
            errorString = f"Variable name \"{name}\" does not exist in particle data"
            errorString = errorString + "\n" + "Following variable names exist: \n"
            for e in self._fields:
                errorString = errorString + "\t" +e + "\n"
            raise NameError(errorString)

        if len(ind) > 1 and n == -1:
            warnings.warn("Multiple entries with the variable name "
                          + name
                          + " exist. Please specify the variable by index.\n"
                          + "E.g. getfield(\"T\",0) returns the first occurance"
                          + " of the variable T"
                          )
        if n==-1:
            n=0
        
        return ind, n
    
    
    def getTimeSeries(self, propertyName, n=-1):
        """
        Retrieves the time values and each particle input property time series
        already sorted based on time values
        USAGE:
            propertyName: name of the property to retrieve the time series
            n: used in case of properties with the same name, check getFieldIndex function
            
        RETURNS:
            timeValues: 1D array with all the time values tracked
            propertyValues: 2D array with the property values for each time step, for 
            each particle tracked
        """
        
        ind, a = self.getfieldIndex(propertyName, n)
        property_index = ind[a]
        timeValues = self._sortedData[:,0,self._time_index]
        sorted_indices = np.argsort(timeValues)
        timeValues = timeValues[sorted_indices]
        
        propertyValues = self._sortedData[:,:,property_index]
        propertyValues = propertyValues[sorted_indices]
        
        return timeValues, propertyValues