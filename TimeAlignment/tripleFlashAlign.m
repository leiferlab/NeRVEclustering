function [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder)
% take flash, yaml, and highRes data for low res fluor, behavior, and high
% res images and makes maps to global time of all movies. Time = 0 will be
% the first frame of the first volume in the hiRes video.

if nargin==0
    mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('FilterSpec',mostRecent);
dataFolder=dataFolder{1};

end


%get data from low yamls(old) or avi metadata text file (new)
if isempty(dir([dataFolder filesep '*.yaml']))
    [~,fluorAll,bfAll]=AviFlashAlign(dataFolder);
else
    [~,fluorAll,bfAll]=YamlFlashAlign(dataFolder);
end

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    %if imsize not specified, will parse from string
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder);
end

hiResData.frameTime=hiResData.frameTime-hiResData.frameTime(1);
bfAll.frameTime=bfAll.frameTime-bfAll.frameTime(1);
fluorAll.frameTime=fluorAll.frameTime-fluorAll.frameTime(1);


hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

%bandaid solution when there is only one flash in each of the recordings.
%If there is one flash in some, and multiple in the other, you're screwed. 

if length(hiResFlashTime)==1 && ...
        length(bfFlashTime)==1 && ...
        length(fluorFlashTime)==1
hiResFlashTime=hiResFlashTime+[0;10];
bfFlashTime=bfFlashTime+[0;10];
fluorFlashTime=fluorFlashTime+[0;10];
end


[~,most]=max([length(hiResFlashTime),length(bfFlashTime),length(fluorFlashTime)]);

%start with track with most flashes detected
switch most
    case 1
        bestFlashTime=hiResFlashTime;
    case 2
        bestFlashTime=bfFlashTime;
    case 3
        bestFlashTime=fluorFlashTime;
end

if isempty(hiResFlashTime) || isempty(bfFlashTime) || isempty(fluorFlashTime)
    error('tripleFlashAlign:noFlash','Flashes not detected in all of the videos')
end


%compare all to the best video with the most flashes
%start with BrightField movie
f_bfTime=timefit(bestFlashTime,bfFlashTime);
bfAll.frameTime=f_bfTime(bfAll.frameTime);


%align with Fluor movie
f_fluorTime=timefit(bestFlashTime,fluorFlashTime);
fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);


f_hiTime=timefit(bestFlashTime,hiResFlashTime);
hiResData.frameTime=f_hiTime(hiResData.frameTime);


%start time is the start of the first volume, make this time zero
startTime=hiResData.frameTime(hiResData.stackIdx==1);
startTime=startTime(1);

hiResData.frameTime=hiResData.frameTime-startTime;
fluorAll.frameTime=fluorAll.frameTime-startTime;
bfAll.frameTime=bfAll.frameTime-startTime;

%recover all the new times for the flashes
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
bfAll.flashTime=bfFlashTime;

fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
fluorAll.flashTime=fluorFlashTime;

hiResFlashTime=hiResData.frameTime(hiResData.flashLoc);
hiResData.flashTime=hiResFlashTime;


function fit_out=timefit(bestFlashTime,sampleFlashTime)
%align with Fluor movie

[~,bf2fluor]=flashTimeAlign2(bestFlashTime,sampleFlashTime);
flashDiff=sampleFlashTime-bestFlashTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);

if length(sampleFlashTime)>1
    fit_out=fit(sampleFlashTime,bestFlashTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
    if fit_out.p1<.1
        fit_out.p1=1;
        
    end
    
else
    fit_out=cfit(fittype('poly1'),1,-sampleFlashTime+bestFlashTime(bf2fluor));
end

