function mutlipleAVIFlash(dataFolder)


d= dir([dataFolder filesep 'LowMagBrain*']);
if length(d)>1
    aviFolder=uipickfiles('filterspec',dataFolder);
    aviFolder=aviFolder{1};
else
aviFolder=[dataFolder filesep d(1).name];
end
camFiles=dir([aviFolder filesep '*.avi']);
camFiles={camFiles.name}';
flashFiles=cellfun(@(x) strrep(x,'.avi','flashTrack.mat'),camFiles,'uniform',0);

if ~exist(flashFiles{1},'file')
    fluorFlash=findFlash(camFiles{1});
end

