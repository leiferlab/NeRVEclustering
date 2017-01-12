function [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize)
% take flash, yaml, and highRes data for low res fluor, behavior, and high
% res images and makes maps to global time of all movies

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
    if nargin<2
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder);
    else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));

    end
    
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
        bestTime=hiResFlashTime;
        bestAll=hiResData;
    case 2
        bestTime=bfFlashTime;
        bestAll=bfAll;
    case 3
        bestTime=fluorFlashTime;
        bestAll=fluorAll;
end

    
%compare all to the best one
%start with BrightField movie
[~,bf2fluor]=flashTimeAlign2(bestTime,bfFlashTime);
flashDiff=bfFlashTime-bestTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);
f_bfTime=fit(bfFlashTime,bestTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_bfTime.p1<.1
    f_bfTime.p1=1;
end
bfAll.frameTime=f_bfTime(bfAll.frameTime);
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
bfAll.flashTime=bfFlashTime;

%align with Fluor movie
[~,bf2fluor]=flashTimeAlign2(bestTime,fluorFlashTime);
flashDiff=fluorFlashTime-bestTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);
if length(fluorFlashTime)>1
f_fluorTime=fit(fluorFlashTime,bestTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_fluorTime.p1<.1
    f_fluorTime.p1=1;
    fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);

end
fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);

else 
     fluorAll.frameTime= fluorAll.frameTime-fluorFlashTime+bestTime(bf2fluor);
end

fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

fluorAll.flashTime=fluorFlashTime;


[~,bf2fluor]=flashTimeAlign2(bestTime,hiResFlashTime);
flashDiff=hiResFlashTime-bestTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);
f_hiTime=fit(hiResFlashTime,bestTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));
if f_hiTime.p1<.1
    f_hiTime.p1=1;
end
hiResData.frameTime=f_hiTime(hiResData.frameTime);
hiResFlashTime=hiResData.frameTime(hiResData.flashLoc);
hiResData.flashTime=hiResFlashTime;
