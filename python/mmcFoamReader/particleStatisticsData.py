import numpy as np
import glob
import warnings
import os
from scipy import stats

class particleStatisticsData:
    
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
    def _numberOfParticlesPerTimeIndex(self,fileList):
        # Open each file and get the line numbers 
        numParticles = np.zeros(len(fileList))
        for file,ind in zip(fileList,range(len(fileList))):
            with open(file, "rb") as f:
                numParticles[ind] = sum(1 for _ in f) - 1

        return numParticles.astype(int)

    def __init__(self,dirPath,fileName):
        """
        Create an object of class particleStatisticsData
        Load the particle statistics into the class and provide functions
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
        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData(
        >>>     "./mmcStatistics/particleSampling/",
        >>>     "probeName")
        Load only a specific time point of a probe
        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData(
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
            # Get the number of particles stored in each file
            self.numParticles = self._numberOfParticlesPerTimeIndex(self._fileList)

            # Total number of particles
            totalNumParticles = np.sum(self.numParticles)

            # Read first entry to get header
            self._setHeader(self._fileList[0])

            # Pre allocate the space for data 
            self._data = np.zeros((totalNumParticles,len(self._fields)))

            currIndex = 0
            for f in self._fileList:
                print("Processing file: ",f)
                data = np.genfromtxt(f,dtype=float,delimiter='\t',skip_header=1)
                if (data.ndim == 1):
                    data = np.reshape(data,(1,len(self._fields)))
                if np.isnan(data[0,-1]):
                    data = data[:,0:-1]
                self._data[currIndex:currIndex+len(data[:,0]),:] = data[:,:].astype(float)
                currIndex = currIndex+len(data[:,0])
        
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

     
    def getNumberOfParticles(self):
        return len(self._data[:,0])

    def getFileList(self):
        """
        Returns the file list read.
        """
        return self._filelist


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
    
    def getfield(self, name, n=-1):
        """
        Get the field data for the given name
        USAGE:
            Returns first occurance of the fieldname
            getfield("nameOfField")
            If multiple occurances with the same name exist the 
            parameter n has to be given to specify which variable should be 
            returned.
            Note, n is zero based, e.g., getField("z",1) returns the 
            conditioning variable not the physical z coordinate 
        """
        ind, n = self.getfieldIndex(name, n)
        return self._data[:,ind[n]]
        
    def conditionalStatistics(self,xName,valuesName,**kwargs):
        """
        Compute conditional statistics data of the variables stored.

        This is a generalization of a histogram function with additional 
        optional parameters. It uses the scipy stats library for the compuation.
        See also the documentation of scipy.stats.binned_statistics() for 
        further information. 


        Parameters
        ----------
        xName : string or (N,) array like
            Name of the particle data field used to condition, e.g. often 
            for mmcFoam the conditioned variable z.
            Can also provide outside data as an array
        valuesName : string or (N,) array like
            Name of the data on which the statistic will be computed. 
            E.g., the temperature 
            Can also provide outside data as an array
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (100 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
        weight: string or (N,) array like, optional
            Weigh the conditioned data by the provided variable name. If an 
            array is provided it has to have the same length as the values
            and x field.

        Returns
        -------
        binCenter : array
            The center point of the bins
        binAv : array
            The binned (conditioned) average data 
        binStdDev : array
            Standard deviation of the binned (conditioned) data in each bin.

        See Also
        --------
        scipy.stats.binned_statistics

        Examples
        --------

        A basic example to calculate the conditional statistics for the 
        temperature T based on the cond. variable z and particle weights wt.

        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData('dir/path/','fileName')
        >>> binCenter, binAv, binStdDev = pData.conditionalStatistics("z","T",weight="wt")
        """

        # Check input
        if not isinstance(xName,str) and not isinstance(xName,type(np.array(10))):
            raise ValueError("xName must be a string or a numpy array")

        if not isinstance(valuesName,str) and not isinstance(valuesName,type(np.array(10))):
            raise ValueError("valuesName must be a string or a numpy array")


        # Get x field
        if isinstance(xName,str):
            x = self.getfield(xName)
        else:
            x = xName
        
        if isinstance(valuesName,str):
            v = self.getfield(valuesName)
        else:
            v = valuesName

        weight = np.ones(len(x))
        bins = min(100,len(x))
        for key, value in kwargs.items():
            if key == "weight":
                # Check type of value
                if isinstance(value,type(np.array(10))):
                    weight=value
                elif isinstance(value,str):
                    weight = self.getfield(value)
                else:
                    raise ValueError("weight must be a numpy array type")
            if key == "bins":
                bins = value
        
        binAv, bin_edges, _ = stats.binned_statistic(x,v*weight,bins=bins,statistic='mean')
      
        binCenter = np.zeros(len(bin_edges)-1)
        for k in range(len(binCenter)):
            binCenter[k] = 0.5*(bin_edges[k+1]+bin_edges[k])
        
        # Conditional mean
        binValxWt, bin_edges, _ = stats.binned_statistic(x,v*weight,bins=bins,statistic='sum')
        binWt, bin_edges, _ = stats.binned_statistic(x,weight,bins=bins,statistic='sum')
        binAv = np.divide(binValxWt,binWt)
        
        # Conditional deviation
        binValSqrxWt, bin_edges, _ = stats.binned_statistic(x,v*v*weight,bins=bins,statistic='sum')
        binStdDev1 = np.divide(binValSqrxWt,binWt)
        binMeanSqr = np.multiply(binAv, binAv)
        binStdDevSqr = np.subtract(binStdDev1,binMeanSqr)
        binStdDev = np.power(binStdDevSqr, 0.5)

        return binCenter,binAv,binStdDev
    
    def conditionalStatisticsRadial(self,axis,valuesName,**kwargs):
        """
        Compute conditional statistics data of the variables stored for 
        a cylindrical domain.

        This is a generalization of a histogram function with additional 
        optional parameters. It uses the scipy stats library for the compuation.
        See also the documentation of scipy.stats.binned_statistics() for 
        further information. 


        Parameters
        ----------
        axis : array like
            Provide the rotation axis of the cylinder as a list or array like 
            type.
        valuesName : string or (N,) array like
            Name of the data on which the statistic will be computed. 
            E.g., the temperature 
            Can also provide outside data as an array
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (100 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
        weight: string or (N,) array like, optional
            Weigh the conditioned data by the provided variable name. If an 
            array is provided it has to have the same length as the values
            and x field.

        Returns
        -------
        binCenter : array
            The center point of the bins
        binAv : array
            The binned (conditioned) average data 
        binStdDev : array
            Standard deviation of the binned (conditioned) data in each bin.

        See Also
        --------
        scipy.stats.binned_statistics

        Examples
        --------

        A basic example to calculate the conditional statistics for the 
        temperature T based on the cond. variable z and particle weights wt
        for a cylinder along the x-axis.

        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData('dir/path/','fileName')
        >>> binCenter, binAv, binStdDev = pData.conditionalStatistics([1,0,0],"T",weight="wt")
        """

        # Check input
        if not isinstance(axis,list) and not isinstance(axis,type(np.array(10))):
            raise ValueError("axis must be a array-like")

        if not isinstance(valuesName,str) and not isinstance(valuesName,type(np.array(10))):
            raise ValueError("valuesName must be a string or a numpy array")


        # create the normal vector to the axis for the radial direction


        # Normalize axis
        axis = axis/np.sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        
        if isinstance(valuesName,str):
            v = self.getfield(valuesName)
        else:
            v = valuesName

        weight = np.ones(len(v))
        bins = min(100,len(v))
        for key, value in kwargs.items():
            if key == "weight":
                # Check type of value
                if isinstance(value,type(np.array(10))):
                    weight=value
                elif isinstance(value,str):
                    weight = self.getfield(value)
                else:
                    raise ValueError("weight must be a numpy array type")
            if key == "bins":
                bins = value
        

        # Calculate the radial position for all particles
        r = np.zeros(len(v))

        # Get x, y, and z coordinates
        x = self.getfield("x",0)
        y = self.getfield("y",0)
        z = self.getfield("z",0)

        k = 0
        for xi,yi,zi in zip(x,y,z):
            # Subtract component in the axial direction
            vec = np.array([xi,yi,zi]) - np.dot(axis,[xi,yi,zi])*axis
            r[k] = np.sqrt(vec.dot(vec))
            k = k+1

        binAv, bin_edges, _ = stats.binned_statistic(r,v*weight,bins=bins,statistic='mean')
      
        binCenter = np.zeros(len(bin_edges)-1)
        for k in range(len(binCenter)):
            binCenter[k] = 0.5*(bin_edges[k+1]+bin_edges[k])
        
        # Conditional mean
        binValxWt, bin_edges, _ = stats.binned_statistic(r,v*weight,bins=bins,statistic='sum')
        binWt, bin_edges, _ = stats.binned_statistic(r,weight,bins=bins,statistic='sum')
        binAv = np.divide(binValxWt,binWt)
        
        # Conditional deviation
        binValSqrxWt, bin_edges, _ = stats.binned_statistic(r,v*v*weight,bins=bins,statistic='sum')
        binStdDev1 = np.divide(binValSqrxWt,binWt)
        binMeanSqr = np.multiply(binAv, binAv)
        binStdDevSqr = np.subtract(binStdDev1,binMeanSqr)
        binStdDev = np.power(binStdDevSqr, 0.5)

        return binCenter,binAv,binStdDev
        

    def reynoldsAveragedConditionalStatistics(self,xName, valuesName, **kwargs):
        """
        Compute reynolds averaged (time averaged) conditional statistics 
        data of the variables stored. (Conditioned on mixture fraction)

        Parameters
        ----------
        xName : string or (N,) array like
            Name of the particle data field used to condition, e.g. often 
            for mmcFoam the conditioned variable z.
            Can also provide outside data as an array
        valuesName : string
            Name of the data on which the statistic will be computed. 
            E.g., the temperature 
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (100 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
        weight: string or (N,) array like, optional
            Weigh the conditioned data by the provided variable name. If an 
            array is provided it has to have the same length as the values
            and x field.

        Returns
        -------
        binCenter : array
            The center point of the bins
        binAv : array
            The binned (conditioned) average data 
        binStdDev : array
            Standard deviation of the binned (conditioned) data in each bin.

        See Also
        --------
        scipy.stats.binned_statistics

        Examples
        --------

        A basic example to calculate the conditional Reynolds averaged statistics
        for the temperature T based on the cond. variable z and particle weights wt.

        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData('dir/path/','fileName')
        >>> binCenter, binAv, binStdDev = pData.reynoldsAveragedConditionalStatistics("T", bins=100)
        """

        # Check input
        if not isinstance(xName,str) and not isinstance(xName,type(np.array(10))):
            raise ValueError("xName must be a string or a numpy array")

        if not isinstance(valuesName,str) and not isinstance(valuesName,type(np.array(10))):
            raise ValueError("valuesName must be a string or a numpy array")


        # Get x field
        if isinstance(xName,str):
            x = self.getfield(xName)
        else:
            x = xName
        
        if isinstance(valuesName,str):
            v = self.getfield(valuesName)
        else:
            v = valuesName

        weight = np.ones(len(x))
        bins = min(100,len(x))
        for key, value in kwargs.items():
            if key == "weight":
                # Check type of value
                if isinstance(value,type(np.array(10))):
                    weight=value
                elif isinstance(value,str):
                    weight = self.getfield(value)
                else:
                    raise ValueError("weight must be a numpy array type")
            if key == "bins":
                bins = value

        # Calculate the time steps
        timeParticles = self.getfield("time")
        timeEdges = np.sort(np.unique(timeParticles))
        timeDeltas = np.zeros(len(timeEdges))

        # Generate the time deltas for time averaging
        
        # First step
        timeDeltas[0] = 0.50 * (timeEdges[0] + timeEdges[1]) - timeEdges[0] 
        
        # Last step
        timeDeltas[-1] = timeEdges[-1] - 0.50 * (timeEdges[-2] + timeEdges[-1])
        
        # Rest of the steps
        for i in range(1, len(timeDeltas)-1):
            timeDeltas[i] = 0.50 * (timeEdges[i] + timeEdges[i+1]) - 0.50 * (timeEdges[i] + timeEdges[i-1])

        # Calculate the bin center and store particles in bins
        _, bin_edges, binIndex = stats.binned_statistic(x,v*weight,bins=bins,statistic='mean')
      
        binCenter = np.zeros(len(bin_edges)-1)
        for k in range(len(binCenter)):
            binCenter[k] = 0.5*(bin_edges[k+1]+bin_edges[k])
        
        # Sort all particles in the binValue
        binValue = np.zeros(len(binCenter))
        # Store the weight to each particle in the bin
        binWeight = np.zeros(len(binCenter))

        # Store the squared weighted value
        binValSqrxWt = np.zeros(len(binCenter))

        deltaTIndex = np.zeros(len(x),dtype=int)


        for i in range(len(x)):
            deltaTIndex[i] = np.where(timeEdges==timeParticles[i])[0][0]
            # Convert to an int for indexing
            binValue[binIndex[i]-1] = binValue[binIndex[i]-1]+(v[i]*weight[i]*timeDeltas[deltaTIndex[i]])
            binWeight[binIndex[i]-1] = binWeight[binIndex[i]-1]+(weight[i]*timeDeltas[deltaTIndex[i]])
            binValSqrxWt[binIndex[i]-1] = binValSqrxWt[binIndex[i]-1]+(v[i]*v[i]*weight[i]*timeDeltas[deltaTIndex[i]])

        binAv = np.divide(binValue,binWeight)
        
        # Conditional deviation
        binStdDev1 = np.divide(binValSqrxWt,binWeight)
        binMeanSqr = np.multiply(binAv, binAv)
        binStdDevSqr = np.subtract(binStdDev1,binMeanSqr)
        binStdDev = np.power(binStdDevSqr, 0.5)

        return binCenter,binAv,binStdDev
    
    
    def reynoldsAveragedConditionalStatisticsRadial(self, axis, valuesName, **kwargs):
        """
        Compute Reynolds averaged conditional statistics data of the variables stored for 
        a cylindrical domain.

        This is a generalization of a histogram function with additional 
        optional parameters. It uses the scipy stats library for the compuation.
        See also the documentation of scipy.stats.binned_statistics() for 
        further information. 


        Parameters
        ----------
        axis : array like
            Provide the rotation axis of the cylinder as a list or array like 
            type.
        valuesName : string or (N,) array like
            Name of the data on which the statistic will be computed. 
            E.g., the temperature 
            Can also provide outside data as an array
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (100 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
        
        Returns
        -------
        binCenter : array
            The center point of the bins
        binAv : array
            The binned (conditioned) average data 
        binStdDev : array
            Standard deviation of the binned (conditioned) data in each bin.

        See Also
        --------
        scipy.stats.binned_statistics

        Examples
        --------

        A basic example to calculate the Reynolds averaged conditional statistics for the 
        temperature T based on the cond. variable z and particle weights wt
        for a cylinder along the x-axis.

        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData('dir/path/','fileName')
        >>> binCenter, binAv, binStdDev = pData.conditionalStatistics([1,0,0],"T",bins=100)
        """
        
                # Check input
        if not isinstance(axis,list) and not isinstance(axis,type(np.array(10))):
            raise ValueError("axis must be a array-like")

        if not isinstance(valuesName,str) and not isinstance(valuesName,type(np.array(10))):
            raise ValueError("valuesName must be a string or a numpy array")


        # create the normal vector to the axis for the radial direction


        # Normalize axis
        axis = axis/np.sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        
        if isinstance(valuesName,str):
            v = self.getfield(valuesName)
        else:
            v = valuesName

        weight = np.ones(len(v))
        bins = min(100,len(v))
        for key, value in kwargs.items():
            if key == "weight":
                # Check type of value
                if isinstance(value,type(np.array(10))):
                    weight=value
                elif isinstance(value,str):
                    weight = self.getfield(value)
                else:
                    raise ValueError("weight must be a numpy array type")
            if key == "bins":
                bins = value
        

        # Calculate the radial position for all particles
        r = np.zeros(len(v))

        # Get x, y, and z coordinates
        x = self.getfield("x",0)
        y = self.getfield("y",0)
        z = self.getfield("z",0)

        k = 0
        for xi,yi,zi in zip(x,y,z):
            # Subtract component in the axial direction
            vec = np.array([xi,yi,zi]) - np.dot(axis,[xi,yi,zi])*axis
            r[k] = np.sqrt(vec.dot(vec))
            k = k+1

        # Generate the time deltas for time averaging
        timeParticles = self.getfield("time")
        timeEdges = np.sort(np.unique(timeParticles))
        timeDeltas = np.zeros(len(timeEdges))

        # First step
        timeDeltas[0] = 0.50 * (timeEdges[0] + timeEdges[1]) - timeEdges[0] 
        
        # Last step
        timeDeltas[-1] = timeEdges[-1] - 0.50 * (timeEdges[-2] + timeEdges[-1])
        
        # Rest of the steps
        for i in range(1, len(timeDeltas)-1):
            timeDeltas[i] = 0.50 * (timeEdges[i] + timeEdges[i+1]) - 0.50 * (timeEdges[i] + timeEdges[i-1])




        # Calculate the bin center and store particles in bins
        _, bin_edges, binIndex = stats.binned_statistic(r,v*weight,bins=bins,statistic='mean')
      
        binCenter = np.zeros(len(bin_edges)-1)
        for k in range(len(binCenter)):
            binCenter[k] = 0.5*(bin_edges[k+1]+bin_edges[k])
        
        # Sort all particles in the binValue
        binValue = np.zeros(len(binCenter))
        # Store the weight to each particle in the bin
        binWeight = np.zeros(len(binCenter))

        # Store the squared weighted value
        binValSqrxWt = np.zeros(len(binCenter))

        deltaTIndex = [timeEdges.index(i) for i in timeParticles]

        for i in range(len(x)):
            binValue[binIndex[i]-1] = binValue[binIndex[i]-1]+(v[i]*weight[i]*timeDeltas[deltaTIndex[i]])
            binWeight[binIndex[i]-1] = binWeight[binIndex[i]-1]+(weight[i]*timeDeltas[deltaTIndex[i]])
            binValSqrxWt[binIndex[i]-1] = binValSqrxWt[binIndex[i]-1]+(v[i]*v[i]*weight[i]*timeDeltas[deltaTIndex[i]])

        binAv = np.divide(binValue,binWeight)
        
        # Conditional deviation
        binStdDev1 = np.divide(binValSqrxWt,binWeight)
        binMeanSqr = np.multiply(binAv, binAv)
        binStdDevSqr = np.subtract(binStdDev1,binMeanSqr)
        binStdDev = np.power(binStdDevSqr, 0.5)
         
        return binCenter, binAv, binStdDev
    
    def calculatePDF(self, valuesName, **kwargs):
        """
        Compute the probability density function of a certain particle 
        variable

        Parameters
        ----------
        valuesName : string or (N,) array like
            Name of the particle data field used to calculate the 
            PDF, e.g. z (mixture fraction), T (temperature), etc...
            Can also provide outside data as an array
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width bins in the
            given range (100 by default).  If `bins` is a sequence, it defines the
            bin edges, including the rightmost edge, allowing for non-uniform bin
            widths.  Values in `x` that are smaller than lowest bin edge are
            assigned to bin number 0, values beyond the highest bin are assigned to
            ``bins[-1]``.  If the bin edges are specified, the number of bins will
            be, (nx = len(bins)-1).
        range : (min, max), optional
            If range provided, bins must be an int, not an array specifying the bin
            edges (must check this, maybe it's the opposite)
        weight: string or (N,) array like, optional
            Weigh the data by the provided variable name. If an 
            array is provided it has to have the same length as the values
            field.

        Returns
        -------
        binsEdges : array, size (bins+1)
            The edges of the resultant PDF bins
        pdfData : array, size (bins)
            The PDF value of each of the bins
        
        See Also
        --------
        numpy.histogram

        Examples
        --------

        A basic example to calculate the PDF for the 
        temperature T considering particle weights wt.

        >>> from mmcFoamReader import particleStatisticsData
        >>> pData = particleStatisticsData('dir/path/','fileName')
        >>> binEdges, binValues = pData.calculatePDF(valuesName="T", bins=100, range=(300, 2000), weight="wt")
        """
        
        # Check input
        if not isinstance(valuesName, str) and not isinstance(valuesName, type(np.array(10))):
            raise ValueError("valuesName must be a string or a numpy array")

        # Get values field
        if isinstance(valuesName,str):
            v = self.getfield(valuesName)
        else:
            v = valuesName

        weight = np.ones(len(v))
        bins = min(100,len(v))
        range = (min(v), max(v))
        for key, value in kwargs.items():
            if key == "weight":
                # Check type of value
                if isinstance(value,type(np.array(10))):
                    weight=value
                elif isinstance(value,str):
                    weight = self.getfield(value)
                else:
                    raise ValueError("weight must be a numpy array type")
            if key == "bins":
                bins = value
            if key == "range":
                range = value
                
        # Calculation of the PDF
        pdfData, binsEdges  = np.histogram(a=v, bins=bins, range=range, weights=weight)
        pdfData = pdfData / sum(pdfData)
        
        return binsEdges, pdfData
