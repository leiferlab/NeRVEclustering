function dataAll=highResTimeTraceAnalysis(imFolder,row,col)
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
%'Frame Number'  'imagesaved' 'DC Offset'  'Image StDev',
labJackData=importdata([imFolder filesep 'LabJackData.txt']);
labJackData=labJackData.data;
%'FuncGen Voltage'  'Z Sensor'  'Photodiode'  'Camera Trigger'

%load images to find flash, flashes have intensity above 10 sigma
datFile=[imFolder filesep 'sCMOS_Frames_U16_1024x1024.dat'];
datFlashRaw=findDatFlash(datFile,10,row,col);
datFlash=datFlashRaw-nanmean(datFlashRaw);
datFlash=datFlash>(nanstd(datFlash)*10);

%% process labJackData
imageWave=normalizeRange(labJackData(:,4));
imageWave=imageWave<.5;
zWave=labJackData(:,2);
%zWave=smooth(zWave,100);

%find time of falling edge
imageSpikes=[0;diff(single(imageWave))>.2]>0;


photoDiode=labJackData(:,3);
photoDiode=photoDiode-smooth(photoDiode,1000);
photoDiode=photoDiode>(std(photoDiode)*10);
photoDiode=[0;diff(single(photoDiode))>.2]>0;
photoDiode=conv(double(photoDiode),ones(10,1),'same');
photoDiode=photoDiode(1:length(labJackData));

imageZ=zWave(imageSpikes);
photoFlash=photoDiode(imageSpikes);


%% process camFrameData

if nnz(photoFlash)~= nnz(datFlash)
    keyboard
    error('Different number of flashes found in image and photoDiode')
end
photoFlashPos=find(photoFlash);
imageIdx=camFrameData(:,1);

imageFlashPos=imageIdx(datFlash);
if length(photoFlashPos)==length(imageFlashPos)
imageIdxOffset=median(imageFlashPos-photoFlashPos);
else
    imageIdxOffset=median(imageFlashPos(1)-photoFlashPos(1));
end


%a given imageIdx now corresponds to both the image z position and the
%image in the imageIdx


dataAll.Z=imageZ;
dataAll.flashLoc=imageFlashPos;
dataAll.imageIdx=(imageIdx);
dataAll.CameraFrameDataOffset=imageIdxOffset;
dataAll.frameTime=find(imageSpikes)/1000;




