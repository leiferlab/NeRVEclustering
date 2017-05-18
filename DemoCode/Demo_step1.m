%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 1 Demo Code: Centerline Fitting
%
%Demo code is made to give a flavor of each section of the analysis
%pipeline. Most of the codes are designed to be run on a large computing
%cluster, but these scripts allow you to run at least part of each section
%locally. 
%
%This script is made for running the centerline fitting code. This must be
%run after "initializeCLworkspace.m" and the file CLworkspace
%

%% Select Folder being analzed,
disp('Select folder containing data to be analyzed');
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% move relevent input data from outputs folder into datafolder for analysis
destination=[dataFolder filesep 'CLworkspace.mat'];
source=[dataFolder filesep 'OutputFiles' filesep 'CLworkspace.mat'];
copyfile(source,destination)

%% run centerline code

%the centerline code is set up to cut the movie into 16 parts, and then
%runs each of those parts forward and backwards, for a total of 32 runs.
%Each of of the 32 runs makes a .mat file in the CLfolder. All of these
%results are compiled later

%can be set to 1 to display some centerlines, runs much slower
display_flag=1; 

%will iterate from 1:32 in current setup. 
run_number=1;
clusterWormCenterline(dataFolder,run_number,display_flag)

%%
% Now, after all of those are done, we have a a folder CL_files full of
% CL.mat files. Move it into the main directory (these are small 
destination_CLfolder=[dataFolder filesep 'CL_files'];
source_CLfolder=[dataFolder filesep 'OutputFiles' filesep 'CL_files'];
if ~exist(destination_CLfolder,'dir')
    movefile(source_CLfolder,destination_CLfolder)
end


%% Compile centerline files, make a final centerline.mat file.
%this requires the CL_files folder located inside the dataFolder dir. The
%outputs of this file is behaviorAnalysis folder with a centerline.mat file
%inside it. 

compileCenterlines(dataFolder)

