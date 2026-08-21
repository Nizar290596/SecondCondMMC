# mmcStatisticsTools
# Functions to read and process the files written by mmcStatistics

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import math



class eulerianStatisticsData:
    """
    Load and processes eulerian statistic files written by 
    mmcFoam with the eulerianStatistics tool.

    Access coordinates with pos() function and values with values()
    To plot data along a line defined by two points 
    p1 = [ax,rad] and p2 = [ax2,rad2] with each an axial and radial
    component use the plotAlongLine(p1,p2) function
    """
    
    # =======================================================================
    # Protected Functions
    def _readInlineList(self,line):
        # First split of the 
        splitLine = line.split("(")
        splitLine = splitLine[1].split(")")
        splitLine = splitLine[0].split()
        val = []
        for e in splitLine:
            val.append(float(e))
        return val

    def _readList(self,f):
        val = []
        while (True):
            line = f.readline()            
            line = line.strip()

            # If it starts with a bracket it is a list
            if line.startswith("(") and line.endswith(")"):
                val.append(self._readInlineList(line))
            elif line.startswith("("):
                continue
            elif line.endswith(")"):
                break
            else:
                val.append(float(line))
        return val
    
    
    # Calculate the magnitude of a vector
    def _magnitude(self,vector):
        return math.sqrt(sum(pow(element, 2) for element in vector))
    
            
    # =======================================================================

    def __init__(self,fname):
        self._fname = fname
        readData = []
        with open(fname) as f:
            while (True):
                # store position of file ptr before readline
                file_pos = f.tell()
                line = f.readline()
                
                # Check for end of file
                if not line:
                    break
                
                line = line.strip()
                # search for the first bracket
                if line.startswith("("):
                    f.seek(file_pos)
                    val = self._readList(f)
                    readData.append(val)
        self._values = readData[0]
        
        # readData[1] currently is a list of list and needs to be converted to an array
        self._pos = np.empty((len(readData[1]),len(readData[1][0])))
        
        rowI = 0
        for row in readData[1]:
            colI = 0
            for e in row:
                self._pos[rowI,colI] = e
                colI = colI + 1
            rowI = rowI + 1
        
        # Generate the interpolator
        self._interp = LinearNDInterpolator(self._pos,self._values)
        
    
    def __str__(self):
        return f"Eulerian data from file: ",_fname
    
    
    # return an x, y data set to plot the data along a line
    # Line is defined with two points
    def plotAlongLine(self,point1,point2,nPoints=100):
        # Make point into a numpy array
        point1 = np.array(point1)
        point2 = np.array(point2)
        
        
        # create points to interpolate to
        t = np.linspace(0.0,1.0,nPoints)
        q = np.empty((len(t),2))
        x = np.linspace(0.0,1.0,nPoints)
        for i in range(len(t)):
            q[i,:] = point1 + t[i]*(point2-point1)
            x[i] = self._magnitude(t[i]*(point2-point1))
        
        # Find the closest point
        v = np.zeros(nPoints)
        for i in range(len(v)):
            v[i] = self._interp(q[i])

        return x,v
    
    def values(self):
        return self._values
    
    def pos(self):
        return self._pos
