import numpy as np
import sys

folder = sys.argv[1]

framesDetails = np.loadtxt(folder+"framesDetails.txt",skiprows=1).T
#framesSync = np.loadtxt(folder+"other-frameSynchronous.dat",skiprows=1).T
#framesAsync = np.loadtxt(folder+"other-frameAsynchronous.dat",skiprows=1).T

frameIndex = framesDetails[1].astype(np.int32)

# Camera frame details -> just synchronous
CameraFrameData = np.array([
    frameIndex,
    frameIndex-frameIndex[0]+1, 
    5.0*np.ones(framesDetails[1].shape[0]),
    np.ones(framesDetails[1].shape[0])
]).T

np.savetxt(folder+"CameraFrameData.txt",CameraFrameData,header="Total Frames\tSaved Frames\tDC Offset\tImage StDev",fmt=['%d','%d','%.6f','%.6f'])
np.savetxt(folder+"LabJackData.txt",np.array([0]))
