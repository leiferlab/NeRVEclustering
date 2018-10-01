import numpy as np
import sys

folder = sys.argv[1]

framesDetails = np.loadtxt(folder+"framesDetails.dat",skiprows=1).T
framesSync = np.loadtxt(folder+"other-frameSynchronous.dat",skiprows=1).T
framesAsync = np.loadtxt(folder+"other-frameAsynchronous.dat",skiprows=1).T

# Camera frame details -> just synchronous
CameraFrameData = np.array([
    framesDetails[0], 
    np.arange(framesDetails[0].shape[0]), 
    5.0*np.ones(framesDetails[0].shape[0]),
    np.ones(framesDetails[0].shape[0])
]).T

np.savetxt(folder+"CameraFrameData.txt",CameraFrameData,header="Total Frames\tSaved Frames\tDC Offset\tImage StDev")
np.savetxt(folder+"LabJackData.txt",np.array([0]))
