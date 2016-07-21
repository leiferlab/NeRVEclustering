
function [initialIm,stackInfo]=timeTraceAnalysis(imFolder)
%  loads excel file with 3 channels representing output camera
% trigger, peizo stage, and flash trigger. searches large image sequence to
% find flash and align flashes and images. initialIm is a max projection of
% ~100 images from the middle section of the stack, stack info is a
% structure array with 1 struct for each z stack, with fields for image
% names, time and names for each image. 
%% load excel file, take outputs and find peak
% if sparsesearch, try to look for flash smartly, otherwise look through
% whole stack.


sparsesearch=0;
customRoi=1;
paraloop=1;
if paraloop;
    paraP=gcp('nocreate');
  if isempty(paraP)
paraP=parpool('local',8,'IdleTimeout', 600);
    end
end


if nargin==0
imFolder=uigetdir();
end


p=dir([imFolder filesep '*.xlsx']);

if isempty(p)
    p=dir([fileparts(imFolder) filesep '*xlsx']);
p=[fileparts(imFolder) filesep p(1).name];
else
p=p(1).name;
p=[imFolder filesep p];

end

imFiles=dir([imFolder filesep '*.tif']);
output=xlsread(p,1,'A:D');

%%
camIdx=2;
zIdx=1;
flashIdx=3;
time=1:length(output);
%find spike for image sync, will find the largest region spanning the
%middle with no spikes in it, ie last spike before middle to first spike
%after
startWave=normalizeRange(output(:,flashIdx));
startWave=(startWave-.5)>.3;
startWave=[0;diff(startWave)];
nFlash=sum(startWave>0);
startTrigger=find(startWave>0 & time'<length(output)/2,1,'last');
endTrigger=find(startWave>0 & time'>length(output)/2,1,'first');
if isempty(endTrigger);
    endTrigger=length(output);
end

startWave=false(size(startWave));
startWave(startTrigger:endTrigger)=true;



%  find image times from the trigger
imageWave=normalizeRange(output(:,camIdx));
imageWave=imageWave<.5;

% find z positions
zWave=output(:,zIdx);
zWave=smooth(zWave,100);

%find time of falling edge
imageSpikes=[0;diff(single(imageWave))>.2];
% imageFall=[diff(single(imageWave))<-.2;0];
% [x,lags]=xcorr(imageFall,imageRise,100);
% shift=lags(x==max(x(lags'>0)) & lags'>0);
% shift=shift(1);f
% imageSpikes=circshift(imageRise,round(shift/2));
% 





%% load images and find flash in images using user defined ROI

stackSize=length(imFiles);
initialIm=(imread([imFolder filesep imFiles(1).name], 'tif'));
if customRoi
progressbar;
for i=round(stackSize/2):round(stackSize/2)+200;
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    initialIm=max(initialIm,temp);
    progressbar((i-round(stackSize/2))/200);

end
imsize=size(initialIm);
fig=imagesc(initialIm);
display('Select area to find flash');
roiFlash=roipoly;
delete(fig)
else
    roiFlash=true(size(initialIm));
end

% fig=imagesc(initialIm);
% display('Select an cell to try to align');
% roiCell1=roipoly;
% delete(fig)
% 
% fig=imagesc(initialIm);
% display('Select another');
% roiCell2=roipoly;
% close all


%% search for the flash
imFlash=zeros(1,stackSize);
%cell1=imFlash;
%cell2=imFlash;
if sparsesearch
progressbar;
for i=1:10:stackSize/10
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    imFlash(i)=mean(temp(roiFlash));
%     cell1(i)=mean(temp(roiCell1));
%     cell2(i)=mean(temp(roiCell2));
    progressbar(i/stackSize*10);
end
progressbar(1);
flashMean=mean(imFlash(imFlash>0));
flashStd=std(imFlash(imFlash>0));


while ~any(imFlash>(flashMean+5*flashStd))
    i=i+10;
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    imFlash(i)=mean(temp(roiFlash));
    flashMean=mean(imFlash(imFlash>0));
flashStd=std(imFlash(imFlash>0));
    
end
i=find(imFlash==max(imFlash));
for k=max(i-100,0):min(i+100,stackSize)
       temp=(imread([imFolder filesep imFiles(k).name], 'tif'));
    imFlash(k)=mean(temp(roiFlash)); 
        progressbar((k-i+100)/200);

    
end
imFlash(imFlash==0)=flashMean;

else
      parfor_progress(stackSize)
    parfor i=1:stackSize
        temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
        imFlash(i)=trimmean(temp(roiFlash),30);
        %     cell1(i)=mean(temp(roiCell1));
        %     cell2(i)=mean(temp(roiCell2));
        %progressbar(i/stackSize);
        parfor_progress;
    end
    parfor_progress(0);
end


%%
try
    

flashWave=imFlash;
flashWave=flashWave-min(flashWave);
flashWave=flashWave>(mean(flashWave)*5);

startFlash=[0,diff(flashWave)>.5];
startFlash=find(startFlash>.5);
startFlashOff=[0,diff(flashWave)<-.5];
startFlashOff=find(startFlashOff>.5);


if length(startFlash)>1
    endFlash=startFlash(find(startFlash>(stackSize/2),1,'first'));
    startFlashIdx=find(startFlash<(stackSize/2),1,'last');
    startFlash=startFlash(startFlashIdx);
    startFlashOff=startFlashOff(startFlashOff>startFlash);
    startFlashOff=startFlashOff(1);
else
    endFlash=stackSize;
end

images=imFiles(startFlashOff:endFlash);
imageZ=zWave(imageSpikes & startWave);
imageZ=imageZ(startFlashOff-startFlash:end);
traceLength=min(length(images),length(imageZ));
imageZ=imageZ(1:traceLength);
images=images(1:traceLength);
zgrad=medfilt1(gradient(zWave,5),5)>0;
zgrad=zgrad(imageSpikes & startWave);
zgrad=zgrad(startFlashOff-startFlash:end);
zgrad=zgrad(1:length(images));
time=time(imageSpikes & startWave);
time=time(startFlashOff-startFlash:end);
time=time(1:length(images));
time=time-min(time);

stackIdx=[0;cumsum(abs(diff(zgrad)))];

for i=1:max(stackIdx);
    imSelect=stackIdx==i;
    stackInfo(i).z=imageZ(imSelect);
    stackInfo(i).fileNames=images(imSelect);
    stackInfo(i).time=time(imSelect);
end


catch ME
save(fullfile(imFolder,'TimeTraceErrorCaught.mat'))
end

if paraloop
    delete(paraP)
end




    



