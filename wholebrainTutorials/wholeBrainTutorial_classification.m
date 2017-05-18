% This is a tutorial script shows how to create neuron registration vectors
% for a single volume with a set of reference files. We'll set it up as a
% demo in which a single volume from a video is used to match the video
% from which it was derived. 

%Requirements:
%   Neuron Registration Vector encoding and Clustering must be completed
%   (step3). The output files PSref.mat and pointStats_info.mat must be
%   present in the folder. 




%% set up the path
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );

ver=version('-release');
ver=str2double(ver(1:4));
if ver<2017
    error('Please use a version of matlab newer than 2017');
end
    

if contains(hostname,'tigressdata')
    cd /tigress/LEIFER/communalCode/3dbrain/
    path(pathdef)
    fprintf('Adding files to path')
else
    disp(['This code is designed to work on tigressdata.'...
        ' You''re not currently on tigressdata so make sure you have the',...
        ' 3dbrain repo in your path!'])
end




%% select the folder Brainscanner folder, this will have all the things we need
% the folder has all the data from the analysis pipeline saved in various
% files. For the time being, all the file names are set so no additional
% inputs are required. 

dataFolder=uipickfiles('Prompt', 'Select the Brain folder', ...
    'FilterSpec','/tigress/LEIFER/PanNeuronal/testing_sets');
dataFolder=dataFolder{1};


%% For our example, we need a pointStats structure from a single volume image. 
%These are produced on mass during the straightening procedure, and are
%compiled into the pointStats mat files in the BrainScanner folder. They
%are created by the clusterWormStraightening program, which does both
%straightening and segmentation. Here, we'll grab a single volume that was
%produced by straightening to use for our example

% search through the file tree and get one of the images
imageFolder=dir([dataFolder filesep 'CLstraigh*']);
imageFolder=[dataFolder filesep imageFolder(1).name];
imageFiles=dir([imageFolder filesep '*.tif']);
target_volume=299;
image_name=[imageFolder filesep imageFiles(target_volume).name];


%load the image into a 3d matrix
V=stackLoad(image_name);

%do segmentation and get pointStats
ps_sample=volume2ps(V);

%% match this pointStats with the reference volumes in the dataFolder, 
% This function runs the non-rigid pointset registration between the
% coordinates in pointStats and the reference pointStats in PSref.mat
% located in the dataFolder. It will take awhile cause it has to do some
% number of registrations. It will then use the registration vectors
% to classify the neurons. 

ps_sample=classifyPS(ps_sample,dataFolder);
%%
% the output pointStats will have 2 new fields, trackIdx, which represent
% the new index for each neuron. The value in trackIdx will match neurons
% in this ps_sample with the neurons in pointStatsNew which should be in
% dataFolder if the analysis pipeline was completed.



