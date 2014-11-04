function dataAll=highResTimeTraceAnalysisTriangle3(imFolder,row,col)
%highResTimeTraceAnalysis takes folder with cameraframedata.txt,
%labjackdata.txt, and sCMOS_Frames_U16_1024x1024.dat and finds alignments
%for timing.

if nargin<2
row=1024;
col=1024;
end
if nargin==0
imFolder=uipickfiles;
imFolder=imFolder{1};

end

%%
camFrameData=importdata([imFolder filesep 'CameraFrameData.txt']);
camFrameData=camFrameData.data;
%'Frame Number'  'DC Offset'  'Image StDev'
imageIdx=camFrameData(:,1);
saveIdx=camFrameData(:,2);

labJackData=importdata([imFolder filesep 'LabJackData.txt']);
labJackData=labJackData.data;
%FuncGen Voltage	Z Sensor	FxnGen Sync	Camera Trigger	Frame Count


%load images to find flash, flashes have intensity above 10 sigma
datFile=[imFolder filesep 'sCMOS_Frames_U16_1024x1024.dat'];
datFlashRaw=findDatFlash(datFile,row,col,10);
datFlash=datFlashRaw-nanmean(datFlashRaw);
datFlash=datFlash>(nanstd(datFlash)*10);

%% process labJackData
imageWave=normalizeRange(labJackData(:,4));
imageWave=imageWave<.5;
%us func gen trigger signal (square) , smooth it out (makes it traingle) and
%use that to sepearte stacks;
zTrigger=labJackData(:,3);
zTrigger=smooth(zTrigger,200);
zWave=labJackData(:,2);

daqSaveFrame=labJackData(:,5);
saveSpikes=diff(daqSaveFrame)>0;

%find time of falling edge
imageSpikes=[0;diff(single(imageWave))>.2]>0;

% for triangles, find up and down
%sepearte rising (1) and falling(0) part of wave
zgrad=medfilt1(gradient(zTrigger,5),3)>0;
zgrad=zgrad(saveSpikes); 
%cumsum the changes in direction, 
stackIdx=[0;cumsum(abs(diff(zgrad)))];
stackAccum=((accumarray(stackIdx+1,ones(size(stackIdx)))));
% get rid of stacks with less than 10 images more than 200 (possibly
% turning off of the zwave);
badIdx=find(stackAccum<10 |stackAccum>200)-1;
sOld=stackIdx;
stackIdx(ismember(sOld,badIdx))=0;
dstackIdx=diff(stackIdx);
dstackIdx(abs(dstackIdx)>1)=0;
stackIdx=[0;cumsum(dstackIdx)];
stackIdx(ismember(sOld,badIdx))=0;

imageZ=zWave(saveSpikes);
spikeFrameIdx=daqSaveFrame(saveSpikes);
timeAll=find(saveSpikes)/1000; % frame time from 1khz 

% [spikeFrameIdx,ib]=unique(spikeFrameIdx);
% timeAll=timeAll(ib);
% imageZ=spikeZ(ib);
% stackIdx=stackIdx(ib);
% imageZ=interp1(spikeFrameIdx,spikeZ,imageIdx,'nearest');
% timeAll=interp1(spikeFrameIdx,timeAll,imageIdx,'nearest');
% stackIdx=interp1(spikeFrameIdx,stackIdx,imageIdx,'nearest');
timeAll(isnan(timeAll))=0;
% camSelect=ismember(saveIdx,1:length(datFlash));
% timeAll=timeAll(camSelect);
% imageZ=imageZ(camSelect);
% imageIdx=imageIdx(camSelect);
% stackIdx=stackIdx(camSelect);

%% process camFrameData

% if nnz(photoFlash)~= nnz(datFlash)
%     keyboard
%     error('Different number of flashes found in image and photoDiode')
% end

imageFlashPos=find(datFlash);



%a given imageIdx now corresponds to both the image z position and the
%image in the imageIdx
imSTD=camFrameData(:,end);
imSTD=imSTD(diff(saveIdx)>0);
imSTD=imSTD-min(imSTD);

imageIdx=imageIdx(diff(saveIdx)>0);
imageIdx=imageIdx-min(imageIdx)+1;
dataAll.Z=imageZ;
dataAll.flashLoc=imageFlashPos;
dataAll.imageIdx=imageIdx;
dataAll.frameTime=timeAll;
dataAll.stackIdx=stackIdx;
dataAll.imSTD=imSTD;
save([imFolder filesep 'hiResData'],'dataAll');



