# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
from matplotlib.pyplot import cm
import numpy as np

# Import libraries for selecting the case folder
from tkinter import filedialog
from tkinter import *
import os

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Select the path to the case
root = Tk()
root.withdraw()
pathToCase = filedialog.askdirectory()

explodeView = True

# Find out how many processors are located in the case folder
nProcs = 0
for file in os.listdir(pathToCase):
    if "processor" in file:
        nProcs = nProcs + 1


subVolumeSetID = input("Select sub-volume id to read in: ")

# Read the sub-volume file
subVolumes = []
with open(pathToCase + "/subVolume-"+str(subVolumeSetID)+".dat","r") as f:
    lines = f.readlines()
    for line in lines:
        temp = np.fromstring(line, dtype=int, sep=' ')
        subVolumes.append(temp)



# Create for each processor the *.foam file to read in paraview
for k in range(nProcs):
    f = open(pathToCase + "/processor"+str(k)+"/proc-"+str(k)+".foam","w")


# create a new 'OpenFOAMReader'
foamReaders = []
for k in range(nProcs):
    foamReaders.append(
        OpenFOAMReader(
            registrationName="proc-"+str(k)+'.foam', 
            FileName=pathToCase + "/processor"+str(k)+'/proc-'+str(k)+'.foam'
        )
    )

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')



# show data in view
foamDisplays=[]
processorRange = np.linspace(0,nProcs-1,nProcs)
for k in processorRange:
    foamDisplays.append(
        Show(foamReaders[int(k)], renderView1, 'UnstructuredGridRepresentation')
    )
    #foamDisplays[-1].SetRepresentationType('Wireframe')


# Make all processors of the same sub-volume in the same color

# Generate range of colors 
colors = cm.rainbow(np.linspace(0, 1, len(subVolumes)))
for subVolume,color in zip(subVolumes,colors):
    for procI in subVolume:
        foamDisplays[procI].AmbientColor = color[0:3]
        foamDisplays[procI].DiffuseColor = color[0:3]



transformations=[]
transformationsDisplay=[]
subVolumeRange = np.linspace(0,len(subVolumes)-1,len(subVolumes))
if explodeView:
    for subVolume,color,i in zip(subVolumes,colors,subVolumeRange):
        # Create a transformation vector 
        r = 0.005
        phi = 2.0*np.pi*float(i)/len(subVolumes)
        theta = np.pi*float(i)/len(subVolumes)
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        for procI in subVolume:
            transformations.append(
                Transform(
                    registrationName='Transform-'+str(procI), 
                    Input=foamReaders[int(procI)]
                    )
            )
            transformations[-1].Transform.Translate = [x,y,z]
            transformationsDisplay.append(
                Show(transformations[-1],renderView1,'UnstructuredGridRepresentation')
            )
            transformationsDisplay[-1].AmbientColor = color[0:3]
            transformationsDisplay[-1].DiffuseColor = color[0:3]
            #Hide(foamDisplays[int(procI)],renderView1)





# trace defaults for the display properties.
#testfoamDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False)

# update the view to ensure updated data information
renderView1.Update()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1017, 784)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 0.0, 0.007692130795259055]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.0010000000474974513]
renderView1.CameraParallelScale = 0.0017320508898368762

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
