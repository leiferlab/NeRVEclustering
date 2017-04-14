function [bf2fluorIdx,fluorAll,bfAll]=AviFlashAlign(dataFolder)


d= dir([dataFolder filesep 'LowMagBrain*']);
if length(d)>1
    aviFolder=uipickfiles('filterspec',dataFolder);
    aviFolder=aviFolder{1};
    
elseif length(d)==1
aviFolder=[dataFolder filesep d(1).name];
elseif isempty(d)
    error('AviFlashAlign:LowMagMissing',...
'ERROR: LowMagBrain folder not found! Ensure LowMagFolder is in BrainScanner Folder')
end

camFiles=dir([aviFolder filesep '*.avi']);
camFiles={camFiles.name}';
camFiles=cellfun(@(x) fullfile(aviFolder,x),camFiles,'uniform',0);

if isempty(camFiles)
    error('AviFlashAlign:AVIMissing',...
'ERROR: avi files in LowMagBrain folder not found! Check that the avi files are present')
end

flashFiles=cellfun(@(x) strrep(x,'.avi','flashTrack.mat'),camFiles,'uniform',0);

if isempty(flashFiles)
    error('AviFlashAlign:FlashMissing',...
'ERROR: Both flash.mat files missing! Have you run flashFinder?')
elseif length(flashFiles)==1
    error('AviFlashAlign:FlashMissing',...
'ERROR: One of the flash.mat files missing! Have you run flashFinder?')
end

if exist(flashFiles{1},'file')
    fluorFlash=load(flashFiles{1});
    fluorFlash=fluorFlash.imFlash;
else
    fluorFlash=findFlash(camFiles{1});
end

if exist(flashFiles{2},'file')
    bfFlash=load(flashFiles{2});
    bfFlash=bfFlash.imFlash;
else
    bfFlash=findFlash(camFiles{2});
end

%%
     
            bfFlash=bfFlash-smooth(bfFlash,200)';
            fluorFlash=fluorFlash-smooth(fluorFlash,200)';
            
bfFlash=bfFlash-min(bfFlash);
%threshold to find peaks
bfFlashloc=find(bfFlash>(mean(bfFlash)+std(bfFlash)*5));
%try to get rid of doubles
bfFlashloc(diff(bfFlashloc)<3)=[];
fluorFlash=fluorFlash-min(fluorFlash);
fluorFlashloc=find(fluorFlash>(mean(fluorFlash)+(std(fluorFlash)*5)));
fluorFlashloc(diff(fluorFlashloc)<3)=[];

camData=importdata([aviFolder filesep 'CamData.txt']);
time=camData.data(:,2);

bf2fluorIdx=1:length(bfFlash);
bfAll.frameTime=time;
fluorAll.frameTime=time;
bfAll.flashTrack=bfFlash;
fluorAll.flashTrack=fluorFlash;
bfAll.flashLoc=union(bfFlashloc,fluorFlashloc);
fluorAll.flashLoc=union(bfFlashloc,fluorFlashloc);
