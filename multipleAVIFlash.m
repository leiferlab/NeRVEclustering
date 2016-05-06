function multipleAVIFlash(dataFolder)


d= dir([dataFolder filesep 'LowMagBrain*']);
if length(d)>1
    aviFolder=uipickfiles('filterspec',dataFolder);
    aviFolder=aviFolder{1};
else
aviFolder=[dataFolder filesep d(1).name];
end

camFiles=dir([aviFolder filesep '*.avi']);
camFiles={camFiles.name}';
for i=1:length(camFiles)
camFileFull=[aviFolder filesep camFiles{i}];
    flashFiles= strrep(camFileFull,'.avi','flashTrack.mat');


if ~exist(flashFiles,'file')
    fluorFlash=findFlash(camFileFull);
end
end
