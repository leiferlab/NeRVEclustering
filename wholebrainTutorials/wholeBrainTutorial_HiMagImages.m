% This is a tutorial script that shows how to load and look at images in
% the .dat file for the whole brain imaging system. The prereq for this
% module is a BrainScanner dataFolder that has had the timing code and
% alignments run (step0 in the wormAnalysis Pipeline)


% Requirements: The first 2 steps of the anlysis pipeline
%       Timing - run submitWormFlashFinder.py
%       alignments - run the alginments_gui on the alignment video.


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

%% closer look at hiResData
% There are many other useful fields in hiResData related to data
% acquisition. Here we'll plot some
%Stack index indentifies each volume in the recording. The value will tell
%the user which volume a given frame belongs to. 
subplot(2,1,1)
plot(hiResData.frameTime,hiResData.stackIdx)
xlim([0,2])
ylim([0,20]);
xlabel('Frame Index')
ylabel('Volume Index');

%The Z value is the voltage reading from the piezo, you'll notice that at
%each change of direction, the stackIdx is incremented. I'll also plot over
%it the image pixel standard dev for each frame. 
subplot(2,1,2)
[ax,~,~]=plotyy(hiResData.frameTime,hiResData.Z,hiResData.frameTime,hiResData.imSTD);
xlim(ax(1),[0,2]);
xlim(ax(2),[0,2]);
ylim(ax(2),[0,50]);

xlabel('Frame Index')
ylabel('Peizo Voltage Reading');
ylabel(ax(2),'Frame pixel sdev');


%% You'll notice on the lower plot that something looks off. 
%You'd expect the Z voltage and the image std to have the same phase, but
%the image std might be a bit offset. This sucks, and is due to software
%delays for processing and saving images. That means that images are offset
%from where we actually think they are. You can calculate this offset using
%cross correlations, or this little program I wrote. 

%take derivative of Zvoltage and heavily smooth, essentially phase shifting
%the signal by -pi/2
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,100);

%normalize the std
image_std=hiResData.imSTD;
image_std=image_std-mean(image_std);
image_std(image_std>100)=0;

%calculate time lag that maximizes the cross correlation between signals
[ZSTDcorrplot,lags]=crosscorr(abs(zWave),image_std,30);
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%lets plot again, but with this offset

%Stack index indentifies each volume in the recording. The value will tell
%the user which volume a given frame belongs to. 
subplot(2,1,1)
plot(hiResData.frameTime,hiResData.stackIdx)
xlim([0,2])
ylim([0,20]);
xlabel('Frame Index')
ylabel('Volume Index');

%The Z value is the voltage reading from the piezo, you'll notice that at
%each change of direction, the stackIdx is incremented. I'll also plot over
%it the image pixel standard dev for each frame. 
subplot(2,1,2)
[ax,~,~]=plotyy(hiResData.frameTime,hiResData.Z,...
    hiResData.frameTime,circshift(hiResData.imSTD,-zOffset));
xlim(ax(1),[0,2]);
xlim(ax(2),[0,2]);
ylim(ax(2),[0,50]);

xlabel('Frame Index')
ylabel('Peizo Voltage Reading');
ylabel(ax(2),'Frame pixel sdev');

display('better right?')

%% Now lets pull some images from the highmag video. This is easy.

%hopefully there are no other .dat files around, get the .dat files name
sCMOSfile=dir([dataFolder filesep '*.dat']);
sCMOSfile=sCMOSfile.name;

%saved in the name is also information about the size of the image, this
%will pull out that image size.
[row,col]=getdatdimensions(sCMOSfile);

%open up the fileID
Fid=fopen([dataFolder filesep sCMOSfile ] );

%lets pull the 99th frame
target_frame=99;
%move the pointer to the image
% each image is row*col pixels, so row*col*2 bytes (we're using uint18),
% move target_frame timees, from the begining (-1)
status=fseek(Fid,2*target_frame*row*col,-1);
%now read the values
fullImage=fread(Fid,row*col,'uint16',0,'l');
%and reshape into an actual image
fullImage=(reshape(fullImage,row,col));
imagesc(fullImage);
title('This is the worm!')

%Thats still not really useful because the two images are stuck together.
%Luckily, we can seperate them using information from the alignments that
%were made

%% loading alignments
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
% green and red image rectangles for cropping, rect1 is red and rect2 is
% green
rect1=alignments.S2AHiRes.rect1;
rect2=alignments.S2AHiRes.rect2;

%we will also need some of these matlab objects for doing a transformation
%to align the green image to the red image
t_concord=alignments.S2AHiRes.t_concord;
Rsegment=alignments.S2AHiRes.Rsegment;

%% Lets get the image again and extract the red and green images
% data concerning the hiMag video alignment is in S2AHires, which has the
% transformation that will take the activity (green) to the segment
% channel(red). It should also have a background image we can use to
% subtract off the signal.

%lets grab the image again
status=fseek(Fid,2*target_frame*row*col,-1);
fullImage=fread(Fid,row*col,'uint16',0,'l');
fullImage=(reshape(fullImage,row,col));
fullImage=fullImage-alignments.background;

%cropping out the red and green image
redImage=fullImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
greenImage=fullImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));

%transforming the green image to match the red image
greenImage=imwarp(greenImage,t_concord,'OutputView',Rsegment);

subplot(2,1,1);
imagesc(redImage)
title('Red image')
subplot(2,1,2);
imagesc(greenImage)
title('Green image, transformed to match the red')

%% Getting an entire volume
% Pulling out an entire volume is very similar, we just need to ask for
% more pixels, and know which pixels we actually want to grab. 
%lets grab the image again

%lets pull out the 99th volume, this will be all the frames that have
starget_stack=99;

target_frames=find(hiResData.stackIdx==starget_stack);

%move the pointer to the first frame
status=fseek(Fid,2*target_frames(1)*row*col,-1);
%now, rather than just row*col pixels, we want
%row*col*length(target_frames)
nPix=row*col*length(target_frames);
fullVolume=fread(Fid,nPix,'uint16',0,'l');
fullVolume=(reshape(fullVolume,row,col,length(target_frames)));
fullVolume=fullVolume-alignments.background;

%cropping out the red and greenf image
redImage=fullVolume((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
greenImage=fullVolume((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);

%transforming the green image to match the red image
greenImage=imwarp(greenImage,t_concord,'OutputView',Rsegment);

%now we can display both of the volumes
h=SliceBrowser(redImage);
h.Name='Red image';
h2=SliceBrowser(greenImage);
h2.Name= 'transformed to match the red';


