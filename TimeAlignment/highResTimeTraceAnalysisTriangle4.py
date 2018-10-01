import numpy as np
import sys
import os

folder = sys.argv[1]

chunkMaxSize = 3000 #frames

filename = folder + "sCMOS_Frames_U16_1024x512.dat"
filesize = os.stat(filename).st_size
f = open(filename)

# Find out how many iterations you have to do
nFrames = filesize // 1024*512*16
nIterations = nFrames // chunkMaxSize
remainingframes = nFrames - nIterations * chunkMaxSize

brightness = np.zeros(nFrames)

for i in nIterations:
    chunk = np.fromfile(f,dtype=np.uint16,count=chunkMaxSize)
    frames = chunk.reshape(chunkMaxSize,1024*512)
    brightness[i*chunkMaxSize:(i+1)*chunkMaxSize] = np.average(frames,axis=1)
    f.seek(chunkMaxSize,1)

chunk = np.fromfile(f,dtype=np.uint16,count=remainingframes)
frames = chunk.reshape(chunkMaxSize,1024*512)
brightness[-remainingframes-1:-1] = np.average(frames,axis=1)

f.close()
plt.plot(brightness)
plt.show()

'''
dataAll = {
    'imageIdx': frameIdx,
    'frameTime': frameIdx*5e-3,
    'flashLoc': flashLoc,
    'stackIdx': volumeIndex,
    'imSTD': frameStd,
    'imAvg': frameAvg,
    'xPos': xPos,
    'yPos': yPos,
}

# Salva questo come mat cavolo di file

# Poi come text file
stringa = "NFrames" + str(volumeIndex[-1])'''
