import numpy as np
import sys
import os
import scipy.io as sio

folder = sys.argv[2]+"/"

# Limit memory usage. Load chunks of data
chunkMaxSize = int(sys.argv[1])//6 # number of frames to read at once

##############
## FLASHFINDER
##############

print("Searching the flashes.")


# Open file containing frames
filename = folder + "sCMOS_Frames_U16_1024x512.dat"
filesize = os.stat(filename).st_size

f = open(filename,'rb')
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
    stdev[i*chunkMaxSize:(i+1)*chunkMaxSize] = np.std(frames.astype(np.float64),axis=1)
    #f.seek(chunkMaxSize*1024*512*2*i) #not needed, np.fromfile moves the pointer

# Last partial chunk
chunk = np.fromfile(f,dtype=np.uint16,count=remainingframes*1024*512)
frames = chunk.reshape((remainingframes,1024*512))
brightness[-remainingframes-1:-1] = np.average(frames,axis=1)
stdev[-remainingframes-1:-1] = np.std(frames.astype(np.float64),axis=1)

# Close file containing frames
f.close()

# Threshold brightness to find flashes
brightnessB = brightness-np.average(brightness)
stdevbrightness = np.std(brightnessB)

# If the flash shows up in two or more consecutive frames, this will list a 
# flash in each of them. 
flashLocRepeated, = np.where(brightnessB>stdevbrightness*10)

# Select only first frame of multiple in which the same flash shows up.
flashLoc =  []
nFlashRep = len(flashLocRepeated)
for nf in np.arange(nFlashRep):
    if nf == 0:
        flashLoc.append(flashLocRepeated[nf])
    else:
        if flashLocRepeated[nf] != flashLocRepeated[nf-1]+1:
            flashLoc.append(flashLocRepeated[nf])
    
flashLoc = np.array(flashLoc)

#################
## VOLUME DETAILS
#################

print("Compiling volume details.")

framesDetails = np.loadtxt(folder+"framesDetails.txt",skiprows=1).T
framesSync = np.loadtxt(folder+"other-frameSynchronous.txt",skiprows=1).T
framesAsync = np.loadtxt(folder+"other-frameAsynchronous.txt",skiprows=1).T
utilities = np.loadtxt(folder+"other-volumeMetadataUtilities.txt",skiprows=1).T

frameIdx = framesDetails[1]
frameTime = framesDetails[0]
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
    'imageIdx': (frameIdx-frameIdx[0]+1).reshape((frameIdx.shape[0],1)),
    'frameTime': (frameTime).reshape((frameTime.shape[0],1)),
    'flashLoc': (flashLoc+1).reshape((flashLoc.shape[0],1)),
    'stackIdx': (volumeIndex+1).reshape((volumeIndex.shape[0],1)),
    'imSTD': stdev.reshape((stdev.shape[0],1)),
    'imAvg': brightness.reshape((brightness.shape[0],1)),
    'xPos': xPos.reshape((xPos.shape[0],1)),
    'yPos': yPos.reshape((yPos.shape[0],1)),
    'Z': Z.reshape((Z.shape[0],1))
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
