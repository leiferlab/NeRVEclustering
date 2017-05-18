% This is a tutorial script that shows how to load and look at images in
% the lowmag .avi files and load the corresponding centerlines. It also shows how
% to do some simple analysis with the centerline.

% Requirements: The first 2 steps of the anlysis pipeline
%       Timing - run submitWormFlashFinder.py
%       Centerlines - submitWormAnalysisCenterline.py


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


%% get the lowMag videos files

%Videos always live in the LowMagBrain folder within the Brainscanner
%folder
d= dir([dataFolder filesep 'LowMagBrain*']);
aviFolder=[dataFolder filesep d(1).name];
%their names are hardcoded to be cam0 and cam1
behaviorMovie=[aviFolder filesep 'cam1.avi'];
fluorMovie=[aviFolder filesep 'cam0.avi'];

%get matlabs video objects, these will have simple stats for the videos,
%such as image size and video length.
behaviorVidObj=VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);

%get the number of frames in the video as frame rate times length
n_frames=behaviorVidObj.FrameRate*behaviorVidObj.Duration;
fprintf('There are %d frames \n', n_frames)

%image sizes are in the Width and Heigh fields
image_size=[behaviorVidObj.Width, behaviorVidObj.Height];
fprintf('The images are %d x %d pixels \n',image_size(1),image_size(2))


%% Now lets pull out an image

target_frame=99;
subplot(1,2,1)
%pull out a specific image for the behavior video
behavior_image=read(behaviorVidObj,target_frame);
imagesc(behavior_image)
title('This is a behavior image')

subplot(1,2,2)
%pull out an image for the fluor video
fluor_image=read(fluorVidObj,target_frame);
imagesc(fluor_image)
title('This is a fluor image')
%% Timing alignment
% The first thing we'll do is load the timing data for all of the movies.
% This is done using the program tripleFlashAlign. The resulting structures
% contain relevent. For the new imaging setup, the fluorecent and behavior
% videos are taken at the same rate, so here bfAll and fluorAll will be
% very close to identical. 

[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);

%Each of these strcutres has a the fields frameTime. These are the time
%vectors for each frame of the low and hiRes videos. The clocks have been
%aligned, so to match a low mag frame with a high mag one, find the same
%time in each of these vectors. 



%% loading Centerlines
%centerlines always live in a behavior folder inside the data folder,

%get behavior folder
behaviorFolder=dir([dataFolder filesep 'Behavior*']);
behaviorFolder=behaviorFolder([behaviorFolder.isdir]);
behaviorFolder=[dataFolder filesep behaviorFolder(1).name];

%in behavior folder, get centerline file
centerlineFile=dir([behaviorFolder filesep 'center*']);
centerlineFile=[behaviorFolder filesep centerlineFile(1).name];

%load centerline file and make variables based on field names
centerline_data=load(centerlineFile);

%centerlinedata will have a number of fields. Each centerline is derived
%from a single frame. 

%% lets look at some of the centerlines and some images
centerline=centerline_data.centerline;

%lets look at image 99
target_frame=99;

behavior_image=read(behaviorVidObj,target_frame);
imagesc(behavior_image)
title('This is a behavior image with a centerline on top')
hold on
%%%%% NOTE %%%%%
% PLOTTING MUST HAPPEN WITH THE SECOND COORDINATE FIRST

%centerlines start at frame 1, so they will have the same index as the
%image. 
plot(centerline(:,2,target_frame),centerline(:,1,target_frame),'r');
hold off

%centerlines have been kind of sad and we're working on imporving the
%system. 

%% Lets look at a few other parts of the centerline data
% the wormcentered coordinate system is used to find the eigenworms.
wormcentered=centerline_data.wormcentered;

imagesc(bfAll.frameTime, [0,1],wormcentered)
set(gca,'ydir','normal')
title('Worm Centered coordinates')
xlabel('Time (s)')
ylabel('Position along the worm')

% You'll notice, if the worm is moving, that there are propagating waves
% along the body of the worm. If the waves are moving from head (1) to
% tail, the animal is moving forward, otherwise it is going backwards.

%% set up the path
cd /tigress/LEIFER/communalCode/3dbrain/
path(pathdef)

%% eigenworms
% the eigenworms are not included with the data, but they live in the
% 3dbrain code repo. 

%load the eigenworms, included are the first 6 eigenworms. 
eigenworms=load('eigenWorms.mat') ;
eigenworms=eigenworms.eigbasis;

% for some reasons, the size of the eigenworms is 101, where our
% centerlines are size 100, if we want to decompose the worm's shape into
% its eigenmodes, we project wormcentered onto this eigenbasis.

%% Luckily, we've already done that

% the projections ofto the first 6 eigenmodes are saved with the
% centerlines
eigenProj=centerline_data.eigenProj;

% eigenProj is a 6 x t times matrix with the projections of each centerline
% onto each of the 6 eigenworms. 

plot(bfAll.frameTime,eigenProj(1,:))
hold on
plot(bfAll.frameTime,eigenProj(2,:))
hold off
xlim([0,10])
xlabel('time(s)')
legend('Projection 1','Projection 2')

%% we can also plot the first 2 modes against eachother,
%for nicely moving worms, this makes a circle, but this doesnt always
%happen
time_range=1:200;
plot(eigenProj(1,time_range),eigenProj(2,time_range));
xlabel('pc1');
ylabel('pc2')


