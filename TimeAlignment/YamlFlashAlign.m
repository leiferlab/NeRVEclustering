function [bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder)
% function takes a folder that has YAML and flashtrack mat files for both a
% fluor image and a behavior image and finds the index mapping between them
% using the positions of the flashes. bf2fluorIdx takes an index in the bf
% image and outputs corresponding index in the fluor series

if nargin==0;
    dataFolder=uipickfiles;
    dataFolder=dataFolder{1};
end

dd=dir([dataFolder filesep '*.mat']);

matName={dd.name}';
yamls=regexpi(matName','yaml');
yamls=cellfun(@(x) ~isempty(x), yamls)';
fluors=regexpi(matName','fluor');
fluors=cellfun(@(x) ~isempty(x), fluors)';
flashes=regexpi(matName','flash');
flashes=cellfun(@(x) ~isempty(x), flashes)';

bits=[flashes,fluors,yamls];
caseNumber=bi2de(fliplr(bits));

for iMat=1:length(matName)
    switch caseNumber(iMat);
        case 4
            bfFlash=load([dataFolder filesep matName{iMat}]);
            bfFlash=bfFlash.imFlash;
        case 6
            fluorFlash=load([dataFolder filesep matName{iMat}]);
            fluorFlash=fluorFlash.imFlash;

        case 3
            fluorYaml=load([dataFolder filesep matName{iMat}]);
            fluorYaml=fluorYaml.mcdf;
        case 1
            bfYaml=load([dataFolder filesep matName{iMat}]);
            bfYaml=bfYaml.mcdf;
    end
end



%% find timing in Yaml files

bfFrameTime=[bfYaml.TimeElapsed]';
bfFrameTime=bfFrameTime-min(bfFrameTime);

fluorFrameTime=[fluorYaml.TimeElapsed]';
fluorFrameTime=fluorFrameTime-min(fluorFrameTime);

%% find flash indecies

bfFlash=bfFlash-min(bfFlash);
bfFlashloc=find(bfFlash>(max(bfFlash)/2));
fluorFlash=fluorFlash-min(fluorFlash);
fluorFlashloc=find(fluorFlash>(max(fluorFlash)/2));

bfFlashTime=bfFrameTime(bfFlashloc);
fluorFlashTime=fluorFrameTime(fluorFlashloc);

bfFlashInterval=diff(bfFlashTime);
fluorFlashInterval=diff(fluorFlashTime);

intervalDif=pdist2(bfFlashInterval,fluorFlashInterval);

[min1,min2]=find(intervalDif==min(intervalDif(:)));

timeDif=bfFlashTime(min1)-fluorFlashTime(min2);

fluorFrameTime=fluorFrameTime+timeDif;


%%
bf2fluorIdx=round(interp1(fluorFrameTime,1:length(fluorFrameTime),bfFrameTime));
fluorAll.frameTime=fluorFrameTime;
fluorAll.flashLoc=fluorFlashloc;
bfAll.frameTime=bfFrameTime;
bfAll.flashLoc=bfFlashloc;

