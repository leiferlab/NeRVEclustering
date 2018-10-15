import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import sys

folder = '/tigress/LEIFER/PanNeuronal/20171023/BrainScanner20171023_162221/'
if os.path.isfile(folder+'hiResData_copy.mat'):
    sys.exit("hiResData.mat has already been corrected for "+folder)
hiResData = sio.loadmat(folder+'hiResData.mat')
sio.savemat(folder+'hiResData_copy.mat',hiResData)

std = np.array(hiResData['dataAll']['imSTD'][0][0][:,0])
stackIdx = np.array(hiResData['dataAll']['stackIdx'][0][0][:,0])

nStack = np.max(stackIdx)
nFrames = len(std)
M = int(round(nFrames/nStack))

deltaStack = stackIdx[1:]-stackIdx[:-1]
stackShift = np.array([i for i, j in enumerate(deltaStack) if j == 1.])+1

#Roll the point of shift between volumes down stdv
for i in np.arange(nStack)[1:]:
    j = stackShift[i]
    dstd = std[j+1]-std[j] #basic gradient descent
    if dstd > 0.:    
        b = np.argmin(std[j-int(M//2):j])
        j = j-int(M//2)+b
    elif dstd < 0.:
        b = np.argmin(std[j:j+int(M//2)])
        j = j+b
    stackShift[i] = j

deltaStackCorrected = np.zeros_like(deltaStack)
deltaStackCorrected[stackShift-1] = 1.0
stackIdxCorrected = np.cumsum(deltaStackCorrected)

hiResData['dataAll']['stackIdx'][0][0][:,0] = np.append(stackIdxCorrected,\
                stackIdxCorrected[-1]).tolist()
sio.savemat(folder+'hiResData.mat',hiResData)

#plt.figure(1,figsize=(15,5))
#plt.plot(std/80.)
#plt.plot(deltaStack)
#plt.plot(stackShift,np.ones(nStack),'ro')
#plt.plot(np.array([i for i, j in enumerate(deltaStack) if j == 1.])+1, 
#            np.ones(nStack),'go')
#plt.plot(deltaStackCorrected,'r-')
#plt.ylim(0,1.1)
#plt.plot(stackIdx)
#plt.plot(stackIdxCorrected)
#plt.xlim(44000,44300)
#plt.xlim(0,1000)
#plt.show()
