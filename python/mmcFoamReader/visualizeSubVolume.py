#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
try:
    from paraview.simple import *
    paraview_module_exist = True
except:
    paraview_module_exist = False
from matplotlib.pyplot import cm
import numpy as np

# Import libraries for selecting the case folder
import os

def visualizeSubVolumeInParaView(pathToCase,**kwargs):
    """
    Visualize the sub-volumes in ParaView

    Parameters
    ----------
    pathToCase : string
        Path to the case folder containing the system/, constant/,
        and processorXX folders
    explodeView: bool (optional)
        Gives an explosion view of the sub-volumes using the transformation 
        function of paraview
    r : scalar (optional)
        Optional parameter to set the radius of the explosion view. 
        By default it is 1.0m

    Examples
    --------
    >>> from mmcFoamReader import visualizeSubVolumeInParaView
    >>> visualizeSubVolumeInParaView('/home/user/OpenFOAM/Case')

    Execute this script with pvpython or load it as a state in ParaView
    """
    explodeView = False

    # Radius for explosion view
    r = 1

    # Process input arguments
    for key,value in kwargs.items():
        if key=="explodeView":
            explodeView = value
        if key=="r":
            r = value


    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    

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
        foamReaders[-1].MeshRegions = ['internalMesh']
        foamReaders[-1].CellArrays = ['']

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
        # First calculate the center point of the complete domain
        centerPointDomain = np.zeros(3)
        centerPointSubVolume = []
        nProcessorDomain  = 0
        for subVolume in subVolumes:
            centerPointSubVolume.append(np.zeros(3))
            for procI in subVolume:
                bounds = foamReaders[int(procI)].GetDataInformation().GetBounds()
                center = np.zeros(3)
                center[0] = 0.5*(bounds[1]+bounds[0])
                center[1] = 0.5*(bounds[2]+bounds[3])
                center[2] = 0.5*(bounds[4]+bounds[5])
                centerPointDomain = centerPointDomain + center
                centerPointSubVolume[-1] = centerPointSubVolume[-1] + center
                nProcessorDomain = nProcessorDomain + 1
            centerPointSubVolume[-1] = centerPointSubVolume[-1]/len(subVolume)
        
        centerPointDomain = centerPointDomain/nProcessorDomain


        # Calculate transformation vector by extending in the direction
        # to the subVolumeCenter

        for subVolume,color,i in zip(subVolumes,colors,subVolumeRange):
            print("Processing subvolume ",i)

            # Calculate vector direction
            dirVec = centerPointSubVolume[int(i)] - centerPointDomain
            
            magSqrDirVec = dirVec[0]**2+dirVec[1]**2+dirVec[2]**2
            if magSqrDirVec == 0:
                dirVec = r*np.ones(3)/np.sqrt(3)
            else:
                # Normalize vector and multiply with r
                dirVec = r*dirVec/(np.sqrt(magSqrDirVec))
            
            for procI in subVolume:
                transformations.append(
                    Transform(
                        registrationName='Transform-'+str(procI), 
                        Input=foamReaders[int(procI)]
                        )
                )
                transformations[-1].Transform.Translate = [dirVec[0],dirVec[1],dirVec[2]]
                transformationsDisplay.append(
                    Show(transformations[-1],renderView1,'UnstructuredGridRepresentation')
                )
                transformationsDisplay[-1].AmbientColor = color[0:3]
                transformationsDisplay[-1].DiffuseColor = color[0:3]
            print("\t Done")





    # trace defaults for the display properties.
    #testfoamDisplay.Representation = 'Surface'

    # reset view to fit data
    #renderView1.ResetCamera(False)

    # update the view to ensure updated data information
    #renderView1.Update()

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================


    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).