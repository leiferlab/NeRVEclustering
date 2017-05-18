%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 5 Demo Code: Signal Extraction
%
%Demo code is made to give a flavor of each section of the analysis
%pipeline. Most of the codes are designed to be run on a large computing
%cluster, but these scripts allow you to run at least part of each section
%locally. 
%
%This script is made for running the Signal extraction after all of the
%analysis pipeline has been completed. It takes the points from the
%straightened coordinates, and then gets the correct pixels from the
%unstraightened red and green channels in order to calculate neural
%signals. The code also does some photobleaching removal. 
%

%% Select Folder being analyzed,
disp('Select folder containing data to be analyzed');
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% move relevent input data from outputs folder into datafolder for analysis

%move the pointStats2 file, which has the results from the registration
%vector encoding, including the preliminary clustering results. 

destination_ps=[dataFolder filesep 'pointStatsNew.mat'];
source_ps=[dataFolder filesep 'OutputFiles' filesep 'pointStatsNew.mat'];
copyfile(source_ps,destination_ps)

destination_behavior=[dataFolder filesep 'BehaviorAnalysis'];
source_behavior=[dataFolder filesep 'OutputFiles' filesep 'BehaviorAnalysis'];
if ~exist(destination,'dir')
    disp('Copying file')
    copyfile(source_behavior,destination_behavior)
else
    disp('File already present');
end

%Error correction requires image intensities, so we also need to tif files
%here. This  is a big folder. You probably want to move this manually
%rather than using matlabs movefile, which can be slow. 
destination_botcheck=[dataFolder filesep 'botCheckFolder'];
source_botcheck=[dataFolder filesep 'OutputFiles' filesep 'botCheckFolder'];
if ~exist(destination_straight,'dir')
    movefile(source_botcheck,destination_botcheck)
end

destination_straight=[dataFolder filesep 'CLstraight'];
source_straight=[dataFolder filesep 'OutputFiles' filesep 'CLstraight'];
if ~exist(destination_straight,'dir')
    movefile(source_straight,destination_straight)
end


%% Doing Signal Extraction
%In practice, this section is looped over the startIdx values as described
%below. This code runs the 

fiducialCropper3(dataFolder)

