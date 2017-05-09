function multipleAVIFlash(dataFolder)
%%% multipleAVIFlash align takes as inpute a BrainScanner folder. It finds
%%% the cam1 and cam0 .avi files int he LowMagBrain folder and finds the
%%% locations of the flashes. 

%find the folder that has the .avi files in it
d= dir([dataFolder filesep 'LowMagBrain*']);
if length(d)>1
    aviFolder=uipickfiles('filterspec',dataFolder);
    aviFolder=aviFolder{1};
else
    aviFolder=[dataFolder filesep d(1).name];
end

%get the avi file names
camFiles=dir([aviFolder filesep '*.avi']);
camFiles={camFiles.name}';
%for each of the avi files, run the findFlash program
for i=1:length(camFiles)
    camFileFull=[aviFolder filesep camFiles{i}];
    flashFiles= strrep(camFileFull,'.avi','flashTrack.mat');

    if ~exist(flashFiles,'file')
        fluorFlash=findFlash(camFileFull);
    end
end
