%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 2 Demo Code: Straightening & Segmentation
%
%Demo code is made to give a flavor of each section of the analysis
%pipeline. Most of the codes are designed to be run on a large computing
%cluster, but these scripts allow you to run at least part of each section
%locally. 
%
%This script is made for running the centerline fitting code. This must be
%run after centerlines
%

%% Select Folder being analyzed,
disp('Select folder containing data to be analyzed');
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% move relevent input data from outputs folder into datafolder for analysis
destination_behavior=[dataFolder filesep 'BehaviorAnalysis'];
source_behavior=[dataFolder filesep 'OutputFiles' filesep 'BehaviorAnalysis'];
if ~exist(destination,'dir')
    disp('Copying file')
    copyfile(source_behavior,destination_behavior)
else
    disp('File already present');
end
%% run straightening code, start
% The straightening code needs to be run on 1 of the volumes first. This is
% done in clusterStraightenStart. This allows subsequent volumes to be
% aligned to the same volume. The output of this program is a
% startWorkspace.mat, which has the information that can be used in the
% next part. There is a copy of this in the output Folder, but as it is
% this doesn't take very long on one core. 

clusterStraightenStart(dataFolder)

%% run straightening code, the rest
% cluster straighten start runs on a range , starting at first_volume,
% going to first_volume+range. It requires a startWorkspace.mat inside the
% dataFolder directory, so it should be run after clusterStraightenStart.
% The output is saved in a CLstraight folder. Currently, tifs and
% pointStats*.mat files are saved, containing each tif . This will be run
% on all volumes in the recording. Currently each volume is being saved as
% a seperate tiff. This is needed later, but can be reprogrammed to save a
% lot of time. 

first_volume=1; %start at begining
range=5;       %lets do 5

clusterWormStraightening(dataFolder,first_volume,range);



%% compile straightening results, prepping them for input into next step
% this takes all of the pointStats*.mat files and compiles all of their
% results and saves a new PointStats.mat file into the dataFolder. This
% requires the CLstraight* Folder in the dataFolder. There are a lot of
% things in the CLstraight folder, you may want to move the folder manually
% from outputFiles to dataFolder. 
msgbox('Move the CLStraight folder from the file output into dataFolder!')

%%
combinePointStatsFiles(dataFolder)


