%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 3 Demo Code: Building Registration Vectors
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
destination_ps=[dataFolder filesep 'PointsStats.mat'];
source_ps=[dataFolder filesep 'OutputFiles' filesep 'PointsStats.mat'];
copyfile(source_ps,destination_ps)

%% Creating registration vectors,
% The program clusterWormTracker builds parts of the registration vectors
% by doing non-rigid point-set registration between sample volumes and
% reference volumes. 

%%%% Inputs
% filePath : link to the complete PointStats file with data from all
% volumes

% startIdx : the index of the the run. If length(pointStats)=N, the total
% number of runs is N*nGroups. For a complete analysis, startIdx will take
% a value 1:N*nGroups for each run

% nGroups : number of runs to each volume. A volume will be compared to
% 150 * nGroups reference volumes

% offset : no longer used
% do group : number of runs to do, will run startIdx: startIdx + doGroups

%removed time limit, no longer needed on della

pointstats_filePath=destination_ps;
startIdx=1;
nGroups=2;
offset=0;
doGroups=1; 


clusterWormTracker(pointstats_filePath,startIdx,nGroups,offset,doGroups)

%% Compile results and cluster each of the registration vectors
% This runs on all of the outputs from the clusterWormTracker, located in
% the folder TrackMatrix. The TrackMatrix folder must be present, along
% with the pointStats.mat file, specified by pointstats_filePath. 

%% move relevent input data from outputs folder into datafolder for analysis
% matlab is very slow at moving files, so you are likely better off moving
% or copying the file in your own file browser. 

destination_trackMatrix=[dataFolder filesep 'TrackMatrix'];
source_trackMatrix=[dataFolder filesep 'OutputFiles' filesep 'TrackMatrix'];
msgbox('Move the TrackMatrix folder from the file output into dataFolder!')

if ~exist(destination_trackMatrix,'dir')
    movefile(source_trackMatrix,destination_trackMatrix)
end

%% also grab the pointStatsFolder that was created after straightening.
destination_ps=[dataFolder filesep 'PointsStats.mat'];
source_ps=[dataFolder filesep 'OutputFiles' filesep 'PointsStats.mat'];
copyfile(source_ps,destination_ps)

%% This step does the clutering of the registration vectors. 
%Both the heiarichical clustering and classification occur here. This is
%computationally intensive and may take up too much memory. 
output_ps=[dataFolder filesep 'PointStats2.mat'];
clusterWormTrackCompiler(source_ps,output_ps)

