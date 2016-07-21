function dataAll=highResTimeTraceAnalysisTriangle(imFolder,row,col)
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
%'FuncGen Voltage'  'Z Sensor'  'Photodiode'  'Camera Trigger' ,'FrameIdx'

%load images to find flash, flashes have intensity above 10 sigma
datFile=[imFolder filesep 'sCMOS_Frames_U16_1024x1024.dat'];
datFlashRaw=findDatFlash(datFile,row,col,10);
datFlash=datFlashRaw-nanmean(datFlashRaw);
datFlash=datFlash>(nanstd(datFlash)*10);

%% process labJackData
imageWave=normalizeRange(labJackData(:,4));
imageWave=imageWave<.5;
zWave=labJackData(:,2);
%zWave=smooth(zWave,100);

daqFrame=labJackData(:,5);

%find time of falling edge
imageSpikes=[0;diff(single(imageWave))>.2]>0;

% for triangles, find up and down
zgrad=medfilt1(gradient(smooth(zWave,10),5),3)>0;
firstRamp=find(zWave>.5,1,'first');
zgrad(1:firstRamp)=0;
zgrad=zgrad(imageSpikes);
stackIdx=[0;cumsum(abs(diff(zgrad)))];



spikeZ=zWave(imageSpikes);
spikeFrameIdx=daqFrame(imageSpikes);
timeAll=find(imageSpikes)/1000; % frame time from 1khz 
[spikeFrameIdx,ib]=unique(spikeFrameIdx);
timeAll=timeAll(ib);
spikeZ=spikeZ(ib);
stackIdx=stackIdx(ib);
imageZ=interp1(spikeFrameIdx,spikeZ,imageIdx,'nearest');
timeAll=interp1(spikeFrameIdx,timeAll,imageIdx,'nearest');
stackIdx=interp1(spikeFrameIdx,stackIdx,imageIdx,'nearest');
timeAll(isnan(timeAll))=0;
camSelect=ismember(saveIdx,1:length(datFlash));
timeAll=timeAll(camSelect);
imageZ=imageZ(camSelect);
imageIdx=imageIdx(camSelect);
stackIdx=stackIdx(camSelect);

%% process camFrameData

% if nnz(photoFlash)~= nnz(datFlash)
%     keyboard
%     error('Different number of flashes found in image and photoDiode')
% end

[ix,ib]=ismember(saveIdx,find(datFlash));

imageFlashPos=imageIdx(ix);



%a given imageIdx now corresponds to both the image z position and the
%image in the imageIdx


dataAll.Z=imageZ;
dataAll.flashLoc=imageFlashPos;
dataAll.imageIdx=(imageIdx);
dataAll.frameTime=timeAll;
dataAll.stackIdx=stackIdx;
save([imFolder filesep 'hiResDataAll'],dataAll);



