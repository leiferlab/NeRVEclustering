%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 4 Demo Code: Error Correction 
%
%Demo code is made to give a flavor of each section of the analysis
%pipeline. Most of the codes are designed to be run on a large computing
%cluster, but these scripts allow you to run at least part of each section
%locally. 
%
%This script is made for running the Error Correction, 
%

%% Select Folder being analyzed,
disp('Select folder containing data to be analyzed');
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% move relevent input data from outputs folder into datafolder for analysis

%move the pointStats2 file, which has the results from the registration
%vector encoding, including the preliminary clustering results. 

destination_ps=[dataFolder filesep 'PointsStats2.mat'];
source_ps=[dataFolder filesep 'OutputFiles' filesep 'PointsStats2.mat'];
copyfile(source_ps,destination_ps)

%Error correction requires image intensities, so we also need to tif files
%here. This  is a big folder. You probably want to move this manually
%rather than using matlabs movefile, which can be slow. 
destination_straight=[dataFolder filesep 'CLstraight'];
source_straight=[dataFolder filesep 'OutputFiles' filesep 'CLstraight'];
if ~exist(destination_straight,'dir')
    movefile(source_straight,destination_straight)
end

%% Doing Error Correction
%In practice, this section is looped over the startIdx values as described
%below. 

%Inputs:
% filepath : path to pointstats2 file made after clusterTrackCompiler
% startIdx : the number of the run. The program is run for each neuron
% cluster detected in the worm. with <groupSize> number of startIdx devoted
% for each neuron
% groupSize: number of groups to run for each neuron. For example, if a
% video is 2000 volumes long, each job will run 2000/groupSize
% comparison. For 150 neurons, the program will need to be run with
% startIdx up to 300 (150 neurons x 2 jobs per neuron) to do every 
% comparison for every neuron. The total computaiton time doesnt change,
% but it just decreases the lenght of each job for more parallelization.


pointstats_filePath=destination_ps;
startIdx=1;
groupSize=2;

clusterBotChecker(pointstats_filePath,startIdx,groupSize)

%% Compile error checking results and do the interpolation 
% This runs on all of the outputs from the clusterBotChecker, located in
% the folder BotChecker. The BotChecker folder, the CLStraight folder, and
% the pointStats2.mat file must be in the dataFolder dir. It saves as an
% output pointStatsNew.mat, which has the refined positions of neurons
% after error correction. This is the final product or the tracking. 

%% move relevent input data from outputs folder into datafolder for analysis
destination_botcheck=[dataFolder filesep 'botCheckFolder'];
source_botcheck=[dataFolder filesep 'OutputFiles' filesep 'botCheckFolder'];
if ~exist(destination_straight,'dir')
    movefile(source_botcheck,destination_botcheck)
end

destination_ps=[dataFolder filesep 'PointsStats2.mat'];
source_ps=[dataFolder filesep 'OutputFiles' filesep 'PointsStats2.mat'];
copyfile(source_ps,destination_ps)

destination_straight=[dataFolder filesep 'CLstraight'];
source_straight=[dataFolder filesep 'OutputFiles' filesep 'CLstraight'];
if ~exist(destination_straight,'dir')
    movefile(source_straight,destination_straight)
end

%%

clusterBotCheckCompiler(dataFolder)


