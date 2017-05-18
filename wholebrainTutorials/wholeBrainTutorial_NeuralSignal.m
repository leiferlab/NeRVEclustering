% This is a tutorial script that shows how to load and look at hi mag
% images in the
% the .dat file for the whole brain imaging system. This tutorial goes over
% the final outputs of the wormAnlaysisPipeline, focusing on the signals
% that are in the heatData.mat file. 

% Requirements: Complete the worm analysis pipeline
%       Timing - run submitWormFlashFinder.py
%       heatData - run submitWormAnalysisPipeline.py through to the end.
%               This file is made by the fiducialCropper, which is the last
%               step of the analysis. 


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

%% First we'll load up the heatData file within the dataFolder

heatFile=[dataFolder filesep 'heatData.mat'];
heatData=load(heatFile)

% As you can see, the heatData is full of different matricies, with
% information for all the neurons and all the times. Lets start with the
% signals.

%% raw data
% the raw data is stored in rRaw and gRaw. These are average pixels for
% each detected neuron without any smoothing or photobleaching correction
% (but with background subtraction for the images).

%the timing vector is in hasPointsTime, the names should be updated after
%Ashley graduates. The time vector is the time for each volume. 
time=heatData.hasPointsTime;

% lets look at one of the neuron's red and green signals
target_neuron=10;
subplot(2,1,1)
plot(time,heatData.rRaw(target_neuron,:),'r','linewidth',3);
xlabel('time')
ylabel('red signal')
title('Raw Red and Green Signals')

subplot(2,1,2)
plot(time,heatData.gRaw(target_neuron,:),'g','linewidth',3)
xlabel('time')
ylabel('green signal')


%% Photobleaching correction is applied to each neuron and each color 
% the signal after photobleaching correction is stored in gPhotoCorr and
% rPhotocorr. 

% lets look at one of the neuron's red and green signals
subplot(2,1,1)
plot(time,heatData.rPhotoCorr(target_neuron,:),'r','linewidth',3);
xlabel('time')
ylabel('red signal')
title('Red and Green Signals after photobleaching correction')

subplot(2,1,2)
plot(time,heatData.gPhotoCorr(target_neuron,:),'g','linewidth',3)
xlabel('time')
ylabel('green signal')

%% smoothed red and green signal
% smothing filters are applied to the red and green signals, along with a
% normalization to show deltaF / F0. F0 is calculated as the lower 20th
% percentile value. 

subplot(2,1,1)
plot(time,heatData.R2(target_neuron,:),'r','linewidth',3);
xlabel('time')
ylabel('red signal')
title('Red and Green Processed Signals after photobleaching correction')

subplot(2,1,2)
plot(time,heatData.G2(target_neuron,:),'g','linewidth',3)
xlabel('time')
ylabel('green signal')


%% Ratio Signal
% the rPhotoCorr and the gPhotoCorr are used to make a ratio sigal. The
% Ratio signal is smoothed and normalized to make an delta R / R0, where R0
% is the lower 20th percentile value. 

plot(time,heatData.Ratio2(target_neuron,:),'r','linewidth',3);
xlabel('time')
ylabel('Ratio signal')
title('Ratiometric neural activity')

%% correlation matrix
% the neurons order is fairly arbitrary, and it is difficult to see
% structure in the neurons. The correlation matrix of the Ratio2 signals
% are shown in the matrix acorr. Each entry is the pearson correlation
% between two neurons. 

imagesc(heatData.acorr)
xlabel('neuron ID')
ylabel('neuron ID')
title('Correlation Matrix (unstructured)')

%% clustered correlation matrix
% the entries of the correlation matrix were heirarchically clustered to
% show structure. This ordering is saved in the field cgIdx.You can
% visualized the ordered correlation matrix as follows:

cgIdx=heatData.cgIdx;
%organize the rows and columns of acorr
imagesc(heatData.acorr(cgIdx,cgIdx))
xlabel('clustered neuron ID')
ylabel('clustered neuron ID')
title('Correlation Matrix (clustered)')

%% clustered neural signals
%cgIdx can also be used to organize the neural signals.

imagesc(time, [], heatData.Ratio2(cgIdx,:))
xlabel('time (s) ')
ylabel( 'clustered neuron ID')
title(' Ratiometric neural signal')



%% Neuron coordinates
% coordinates of the neurons in one of the straightened coordinate systems
% is saved in XYZcoord. The coordinates are in the units of pixels in the
% high mag camera
XYZcoord=heatData.XYZcoord;
n_neurons=size(XYZcoord,1);

%we'll show an example scatter plot with each neuron shown by a sphere. 
scatter3sph(XYZcoord(:,1),XYZcoord(:,2),XYZcoord(:,3),'size',5);

%we'll also label all of the balls with their number
text(XYZcoord(:,1),XYZcoord(:,2),XYZcoord(:,3),...
     cellstr(num2str((1:n_neurons)'))) %make cellstr of numbers
    
axis equal

%% behavior

%a summary of behaviors are provided at the same timing as the neural
%recordings. This data is stored in the behavior field

behavior= heatData.behavior;

% the data has several behavior metrics, the simplest one is an ethogram of
% the behavior. 


%the ethogram has -1 for backward, 0, for pause, 1 for forward, 2 for turn

% define ethogramColormap
forward_color=[0 1 0];% green
back_color=[1 0 0];%red
turn_color=[0 0 1];%blue
pause_color=[255 217 50]/256; %yellowish
ethocolormap=[back_color;pause_color;forward_color;turn_color];

% to view the ethogram with colors
imagesc(time,[],behavior.ethogram');
caxis([-1,2])


%% other behavior metrics

%the x and y position are the coordinates of the center of mass of the worm
%in the coordinate system of the plate, in units of mm. 
subplot(2,1,1)
plot(time,behavior.x_pos)
ylabel('x position (mm)')

subplot(2,1,2)
plot(time,behavior.y_pos)
ylabel('y position (mm)')
xlabel('time (s)')
