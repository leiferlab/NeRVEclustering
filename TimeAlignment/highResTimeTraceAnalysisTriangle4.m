function dataAll=highResTimeTraceAnalysisTriangle4(imFolder)
%highResTimeTraceAnalysis takes folder with cameraframedata.txt,
%labjackdata.txt, and sCMOS_Frames_U16_<rows>x<cols.dat and finds alignments
%for timing.
%modified 20160626 to work with dat files with different sizes

%if no input, select folder using uipickfiles
if nargin==0
imFolder=uipickfiles;
imFolder=imFolder{1};

end

%% load text files
camFrameData=importdata([imFolder filesep 'CameraFrameData.txt']);
camFrameData=camFrameData.data;
% camFrameData is text file with columns:
%'Frame Number'  'DC Offset'  'Image StDev'
saveIdx=camFrameData(:,2);


labJackData=importdata([imFolder filesep 'LabJackData.txt']);
labJackData=labJackData.data;
% LabJackData is text file with columns:
%FuncGen Voltage ZSensor FxnGenSync CameraTrigger savedFrameCount FrameCount {x y}

%%
%load images to find flash, flashes have intensity above 10 sigma
datFile=dir([imFolder filesep 'sCMOS_Frames_U16_*.dat']);
datFile=[imFolder filesep datFile.name];

%get size of dat file if image dimensions are not in input.
[row,col]=getdatdimensions(datFile);
%%
%make time trace of intensity to find flash
datFlashRaw=findDatFlash(datFile,row,col,10);
%threshold to find flash.
datFlash=datFlashRaw-nanmean(datFlashRaw);
datFlash=datFlash>(nanstd(datFlash)*10);
datFlash=find(datFlash);

%try to fix fast doubles
datFlash(diff(datFlash)<3)=[];
%% process labJackData
imageWave=normalizeRange(labJackData(:,4));
imageWave=imageWave<.5;
%us func gen trigger signal (square) , smooth it out (makes it traingle) and
%use that to sepearte stacks;
zTrigger=labJackData(:,3);
zTrigger2=normalizeRange(zTrigger)>.5;
smoothKernal=median(diff(find(diff(zTrigger2)==1)))/2;
zTrigger=smooth(zTrigger,smoothKernal);
zWave=labJackData(:,2);

daqSaveFrame=labJackData(:,5);
saveSpikes=diff(daqSaveFrame)>0;

%find time of falling edge
imageSpikes=[0;diff(single(imageWave))>.2]>0;

% for triangles, find up and down
%sepearte rising (1) and falling(0) part of wave
%zTrigger=normalizeRange(zTrigger);
%zTrigger=double(zTrigger>.5);
zgrad=zTrigger; 
zgrad=(gradient(zgrad));
%zgrad=abs([0 ;(diff(zgrad))]);
zgrad=(zgrad/std(zgrad))>.1;
zgrad=[0; diff(zgrad)];
zgrad=abs(zgrad)>.1;
%zgrad=zgrad>max(zgrad/2);
%cumsum the changes in direction, 
stackIdx=[0;cumsum(zgrad)];
stackIdx=stackIdx(saveSpikes);
stackAccum=((accumarray(stackIdx+1,ones(size(stackIdx)))));
% get rid of stacks with less than 10 images more than 200 (possibly
% turning off of the zwave);
badIdx=find(stackAccum<mean(stackAccum*.5) |stackAccum>200)-1;
sOld=stackIdx;
stackIdx(ismember(sOld,badIdx))=0;
dstackIdx=diff(stackIdx);
%dstackIdx((dstackIdx)>1)=0;
stackIdx=[0;cumsum(dstackIdx>0)];
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

imageIdx=camFrameData(:,1);



%a given imageIdx now corresponds to both the image z position and the
%image in the imageIdx
imSTD=camFrameData(:,end);
imSTD=imSTD(diff(saveIdx)>0);
imSTD=imSTD-min(imSTD);

imageIdx=imageIdx(diff(saveIdx)>0);
imageIdx=imageIdx-min(imageIdx)+1;
dataAll.Z=imageZ;
dataAll.flashLoc=datFlash;
dataAll.imageIdx=imageIdx;
dataAll.frameTime=timeAll;
dataAll.stackIdx=stackIdx;
dataAll.imSTD=imSTD;
dataAll.imAvg=datFlashRaw;
if size(labJackData,2)>7
    xPos=labJackData(:,7);
    yPos=labJackData(:,8);
    dataAll.xPos=xPos(saveSpikes);
    dataAll.yPos=yPos(saveSpikes);
else
    dataAll.xPos=[];
    dataAll.yPos=[];
    
end


save([imFolder filesep 'hiResData'],'dataAll');

%write a little text file so python can read in number of stackss
Fid=fopen([imFolder filesep 'submissionParameters.txt'],'w');
nStacks=max(stackIdx);
fprintf(Fid,'%s %d','NFrames',nStacks);
fclose(Fid);


