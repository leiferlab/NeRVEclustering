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

if sum(yamls)<2 ||sum(flashes)<2
    dyamls=dir([dataFolder filesep '*.yaml']);
    if isempty(dyamls)
        if length( yamls)<2
            display('Select Yaml Files');
            yamlFiles=uipickfiles('filterspec', dataFolder);
        end
    else
        yamlFiles={dyamls.name};
        yamlFiles=cellfun(@(x) [dataFolder filesep x],yamlFiles,'uniform',0);
    end
    davis=dir([dataFolder filesep '*.avi']);
    davis={davis.name};
    davis=davis(cellfun(@(x) isempty(strfind(x, 'HUDS')),davis));
    if isempty(davis)
        if length(flashes)<2
            display('Select avi files');
            aviFiles=uipickfiles('filterspec', dataFolder);
        end
    else
        aviFiles=davis;
        aviFiles=cellfun(@(x) [dataFolder filesep x],aviFiles,'uniform',0);
        
    end
    
    findFlash(aviFiles);
    makeMatFileFromYaml(yamlFiles);
    
    dd=dir([dataFolder filesep '*.mat']);
    matName={dd.name}';
    yamls=regexpi(matName','yaml');
    yamls=cellfun(@(x) ~isempty(x), yamls)';
    fluors=regexpi(matName','fluor');
    fluors=cellfun(@(x) ~isempty(x), fluors)';
    flashes=regexpi(matName','flash');
    flashes=cellfun(@(x) ~isempty(x), flashes)';
    
end



bits=[flashes,fluors,yamls];
caseNumber=bi2de(fliplr(bits));
%%
for iMat=1:length(matName)
    switch caseNumber(iMat);
        case 4
            
            bfFlash=load([dataFolder filesep matName{iMat}]);
            bfFlash=bfFlash.imFlash;
            bfFlash=bfFlash-smooth(bfFlash,200)';
        case 6
            fluorFlash=load([dataFolder filesep matName{iMat}]);
            fluorFlash=fluorFlash.imFlash;
            fluorFlash=fluorFlash-smooth(fluorFlash,200)';
            
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
bfFlashloc=find(bfFlash>(mean(bfFlash)+std(bfFlash)*5));
fluorFlash=fluorFlash-min(fluorFlash);
fluorFlashloc=find(fluorFlash>(mean(fluorFlash)+(std(fluorFlash)*5)));

bfFlashTime=bfFrameTime(bfFlashloc);
fluorFlashTime=fluorFrameTime(fluorFlashloc);

bfFlashInterval=diff(bfFlashTime);
fluorFlashInterval=diff(fluorFlashTime);
if ~isempty(fluorFlashInterval) && ~isempty(bfFlashInterval)
    intervalDif=pdist2(bfFlashInterval,fluorFlashInterval);
    
    [min1,min2]=find(intervalDif==min(intervalDif(:)));
    timeDif=bfFlashTime(min1)-fluorFlashTime(min2);
    
else
    timeDif=bfFlashTime(1)-fluorFlashTime(1);
end
fluorFrameTime=fluorFrameTime+timeDif;


%%
bf2fluorIdx=round(interp1(fluorFrameTime,1:length(fluorFrameTime),bfFrameTime));
fluorAll.frameTime=fluorFrameTime;
fluorAll.flashLoc=fluorFlashloc;
fluorAll.flashTrack=fluorFlash;
bfAll.frameTime=bfFrameTime;

bfAll.flashLoc=bfFlashloc;
bfAll.flashTrack=bfFlash;

