
function [initialIm,stackInfo]=timeTraceAnalysis_screwedFlash(imFolder)
%  loads excel file with 3 channels representing output camera
% trigger, peizo stage, and flash trigger. searches large image sequence to
% find flash and align flashes and images. initialIm is a max projection of
% ~100 images from the middle section of the stack, stack info is a
% structure array with 1 struct for each z stack, with fields for image
% names, time and names for each image.


sparsesearch=1;

%% load excel file, take outputs and find peak
if nargin==0
    imFolder=uigetdir();
end

p=dir([imFolder filesep '*.xlsx']);
p=p(1).name;
imFiles=dir([imFolder filesep '*.tif']);
output=xlsread([imFolder filesep p],1,'A:C');

%%
camIdx=1;
zIdx=2;
flashIdx=3;

time=1:length(output);
%find spike for image sync
startWave=normalizeRange(output(:,flashIdx));
startWave=startWave==max(startWave);
startWaveSmooth=normalizeRange(smooth(startWave,500));
startWave=startWave & startWaveSmooth>.5;

startWave=cumsum(startWave>.2)>1;



%  find image times from the trigger
imageWave=normalizeRange(output(:,camIdx));
imageWave=imageWave<.5;

% find z positions
zWave=output(:,zIdx);
zWave=smooth(zWave,100);

%find time of images, take middle of rising and falling edge
imageSpikes=[0;diff(single(imageWave))>.2];
% imageFall=[diff(single(imageWave))<-.2;0];
% [x,lags]=xcorr(imageFall,imageRise,100);
% shift=lags(x==max(x(lags'>0)) & lags'>0);
% shift=shift(1);
% imageSpikes=circshift(imageRise,round(shift/2));
%





%% load images and find flash in images using user defined ROI

stackSize=length(imFiles);
initialIm=(imread([imFolder filesep imFiles(1).name], 'tif'));
progressbar;
for i=round(stackSize/2):round(stackSize/2)+400;
    temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
    initialIm=max(initialIm,temp);
    progressbar((i-round(stackSize/2))/400);
    
end
progressbar(1);
imsize=size(initialIm);
fig=imagesc(initialIm);
display('Select area to find flash');
roiFlash=roipoly;
delete(fig)

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

if sparsesearch
    %cell1=imFlash;
    %cell2=imFlash;
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
    progressbar
    for i=1:1:stackSize
        temp=(imread([imFolder filesep imFiles(i).name], 'tif'));
        imFlash(i)=mean(temp(roiFlash));
        %     cell1(i)=mean(temp(roiCell1));
        %     cell2(i)=mean(temp(roiCell2));
        progressbar(i/stackSize);
    end
    
end


%%

flashWave=imFlash;
flashWave=flashWave-min(flashWave);
flashWave=flashWave>(mean(flashWave)*5);
startFlash=[0,diff(flashWave)>.5];
startFlash=find(startFlash>.5);
startFlashOff=[0,diff(flashWave)<-.5];
startFlashOff=find(startFlashOff>.5);


if length(startFlash)>1
    endFlash=startFlash(end);
    startFlash=startFlash(1);
else
    endFlash=stackSize;
end

images=imFiles(startFlashOff:endFlash);
imageZ=zWave(imageSpikes & startWave);
imageZ=imageZ(startFlashOff-startFlash:end);
imageZ=imageZ(1:length(images));
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















