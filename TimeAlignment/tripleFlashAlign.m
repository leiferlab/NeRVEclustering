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
[~,bf2fluor]=flashTimeAlign2(bestFlashTime,bfFlashTime);
flashDiff=bfFlashTime-bestFlashTime(bf2fluor);

flashDiff=flashDiff-min(flashDiff);
f_bfTime=fit(bfFlashTime,bestFlashTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_bfTime.p1<.1
    f_bfTime.p1=1;
end
bfAll.frameTime=f_bfTime(bfAll.frameTime);


%align with Fluor movie
[~,bf2fluor]=flashTimeAlign2(bestFlashTime,fluorFlashTime);
flashDiff=fluorFlashTime-bestFlashTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);
if length(fluorFlashTime)>1
f_fluorTime=fit(fluorFlashTime,bestFlashTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_fluorTime.p1<.1
    f_fluorTime.p1=1;
    fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);

end
fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);

else 
     fluorAll.frameTime= fluorAll.frameTime-fluorFlashTime+bestFlashTime(bf2fluor);
end




[~,bf2fluor]=flashTimeAlign2(bestFlashTime,hiResFlashTime);
flashDiff=hiResFlashTime-bestFlashTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);
f_hiTime=fit(hiResFlashTime,bestFlashTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_hiTime.p1<.1
    f_hiTime.p1=1;
end
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
