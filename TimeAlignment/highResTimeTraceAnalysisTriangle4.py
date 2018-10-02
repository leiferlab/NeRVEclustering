import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import scipy.io as sio

folder = sys.argv[2]

# Limit memory usage
chunkMaxSize = int(sys.argv[1]) # number of frames to read at once

##############
## FLASHFINDER
##############

print("Searching the flashes.")

# Open file containing frames
filename = folder + "sCMOS_Frames_U16_1024x512.dat"
filesize = os.stat(filename).st_size
f = open(filename)

# Find out how many iterations you have to do
nFrames = filesize // (1024*512*2)
nIterations = nFrames // chunkMaxSize
remainingframes = nFrames - nIterations * chunkMaxSize
brightness = np.zeros(nFrames)
stdev = np.zeros(nFrames)

# Full iterations
for i in np.arange(nIterations):
    chunk = np.fromfile(f,dtype=np.uint16,count=chunkMaxSize*1024*512)
    frames = chunk.reshape((chunkMaxSize,1024*512))
    brightness[i*chunkMaxSize:(i+1)*chunkMaxSize] = np.average(frames,axis=1)
    stdev[i*chunkMaxSize:(i+1)*chunkMaxSize] = np.std(frames, axis=1)
    f.seek(chunkMaxSize*1024*512*2)

# Last partial chunk
chunk = np.fromfile(f,dtype=np.uint16,count=remainingframes*1024*512)
frames = chunk.reshape((remainingframes,1024*512))
brightness[-remainingframes-1:-1] = np.average(frames,axis=1)

# Close file containing frames
f.close()

# Threshold brightness to find flashes
brightnessB = brightness-np.average(brightness)
stdevbrightness = np.std(brightnessB)

flashLoc, = np.where(brightnessB>stdevbrightness*10)

#################
## VOLUME DETAILS
#################

print("Compiling volume details.")

framesDetails = np.loadtxt(folder+"framesDetails.dat",skiprows=1).T
framesSync = np.loadtxt(folder+"other-frameSynchronous.dat",skiprows=1).T
framesAsync = np.loadtxt(folder+"other-frameAsynchronous.dat",skiprows=1).T
utilities = np.loadtxt(folder+"other-volumeMetadataUtilities.dat",skiprows=1).T

frameIdx = framesDetails[1]
volumeIndex = np.zeros_like(framesDetails[1],dtype=np.int32)
xPos = np.zeros_like(framesDetails[1])
yPos = np.zeros_like(framesDetails[1])
Z = framesSync[1]

vindold = -1
vsignold = 0
xpos = 0.0
ypos = 0.0
for i in np.arange(len(framesDetails[1])):
    frameindex = framesDetails[1][i]
    
    # Assign a volume index to each frame
    try:
        j1 = np.where(framesSync[0]==frameindex)[0][0]
        
        if framesSync[2,j1] != vsignold:
            vindold += 1
        vsignold = framesSync[2,j1] 
        volumeIndex[i] = vindold
    except:
        volumeIndex[i] = vindold
    
    # Assign the xpos and ypos from the ludl stage to each frame
    try:
        j2 = np.where(framesAsync[0]==frameindex)[0][0]
        xpos = framesAsync[1][j2]
        ypos = framesAsync[2][j2]
    except:
        pass
         
    xPos[i] = xpos
    yPos[i] = ypos
        
Mat = {}
dataAll = {
    'imageIdx': frameIdx-frameIdx[0]+1,
    'frameTime': (frameIdx-frameIdx[0]+1)*5e-3,
    'flashLoc': flashLoc+1,
    'stackIdx': volumeIndex+1,
    'imSTD': stdev,
    'imAvg': brightness,
    'xPos': xPos,
    'yPos': yPos,
    'Z': Z
}
Mat['dataAll'] = dataAll

sio.savemat(folder + "hiResData.mat",Mat)

# Salva questo come mat cavolo di file


##########################
## OLD CameraFrameData.txt
##########################

print("Printing files.")

# Camera frame details -> just synchronous
CameraFrameData = np.array([
    framesDetails[0], 
    np.arange(framesDetails[0].shape[0]), 
    5.0*np.ones(framesDetails[0].shape[0]),
    np.ones(framesDetails[0].shape[0])
]).T

np.savetxt(folder+"CameraFrameData.txt",CameraFrameData,
            header="Total Frames\tSaved Frames\tDC Offset\tImage StDev")


##########################################
## FILE TELLING HOW MANY VOLUMES THERE ARE
##########################################

stringa = "NFrames " + str(volumeIndex[-1])

f = open(folder + "submissionParameters.txt", "w")
f.write(stringa)
f.close()
