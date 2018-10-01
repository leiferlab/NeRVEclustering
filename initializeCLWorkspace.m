
function initializeCLWorkspace(dataFolder)
% This script is used for initializing parameters and data for centerline
% detection. This script can be run cell by cell, with a few or the cells
% running for or parfor loops over every frame of the video, and some cells
% require user inpurts to initialize centerlines or crop out image regions.
% Prior to running this script, you may run wormCL_tip_clicker.m for a GUI
% that allows the user to click centerline tips.



%% parameters
%% Initialize fitting parameters for centerline, 

cline_para.refIdx=8;
cline_para.tipRegion=45;
cline_para.endRepulsion=.5;
cline_para.repulsionD=20;
cline_para.heat=3;
cline_para.CLalpha=1;
cline_para.CLbeta=40;
cline_para.gamma=40;  
cline_para.kappa=60;
cline_para.endkappa=5;
cline_para.gradient_force=2;
cline_para.showFlag=0;
cline_para.iterations=400;
cline_para.stretching_force_factor=[.3 .3];
cline_para.refSpring=.01;
cline_para.stretch_ends_flag=1;
cline_para.refL=5.5;
cline_para.memForce=.01;


% these are from autoCL and tend to work better. maybe slower
%cline_para.memForce=0;   
%cline_para.stretching_force_factor=[1 0.3];
%cline_para.gradient_force=1;
%cline_para.refSpring=0;
    



%% Select datafolder for analysis
if nargin==0
mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('FilterSpec',mostRecent,...
    'Prompt', 'Select the BrainScanner Folder or the LowMag Folder');
dataFolder=dataFolder{1};
end
setappdata(0,'mostRecent',fileparts(dataFolder));

display(['Data Folder: ', dataFolder]);


aviFiles=dir([dataFolder filesep '*.avi']);
aviFiles={aviFiles.name}';
HUDFiles=aviFiles(cellfun(@(x) ~isempty(strfind(x,'HUDS')),aviFiles));
oldFlag=length(HUDFiles);

if ~oldFlag
    
    %% get lowmag folder
    if strfind(dataFolder,'LowMag')
        low_mag_folder=dataFolder;
    else
        low_mag_folder=dir([dataFolder filesep 'LowMag*']);
        if isempty(low_mag_folder)
            
            error(...
                'LowMag folder is missing! ensure the Low mag folder is in the BrainScanner Folder')
        end
        low_mag_folder=[dataFolder filesep low_mag_folder(1).name];
    end
    %setup paths to movies
    fluormovie=[low_mag_folder filesep 'cam0.avi'];
    behaviormovie=[low_mag_folder filesep 'cam1.avi'];
    %get movie length
    bf2fluor_lookup=[];
else
    low_mag_folder=dataFolder;

    aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviormovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluormovie=[dataFolder filesep aviFiles{aviFluorIdx}];
    %% set up timing alignments and lookups
    
    %get timing sync for old data movies, folders were hard saved at
    %1200x600
    [bfAll,fluorAll,~]=tripleFlashAlign(dataFolder);
    fluorIdxList=1:length(fluorAll.frameTime);
    bf2fluor_lookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');
    
end


%% set up low magvideos, we've changed the way we save data, the older version

%initialize video objects
behavior_vidobj = VideoReader(behaviormovie);
fluor_vidobj= VideoReader(fluormovie);
%get move length
nframes=round(behavior_vidobj.Duration*behavior_vidobj.FrameRate);
nframes_fluor=round(fluor_vidobj.Duration*fluor_vidobj.FrameRate);
bf_imsize=[behavior_vidobj.Height,behavior_vidobj.Width];

if ~oldFlag && nframes~=nframes_fluor
    error('Videos should be of the same length, check if one of the videos is corrupt')
end

%% initialize centerline points

%cut up the entire video into nCells chunks and then initialize nCells+1
%centerlines for fitting each chunk forwards and backwards.
nCells=16;
nSteps=ceil(nframes/nCells); %number of frames in each chunk (except last)
bfCell_i=cell(nCells,1);
clStartI=cell(nCells+1,1);
display(['Click the ~5-10 points of the centerline starting at the head. When '...
    'you are done, double click the last point.']);
%
sample_images=zeros(bf_imsize(1),bf_imsize(1),nCells+1);
for ichunk=1:nCells+1
    %get bounds for each chunk
    lowframe=min((ichunk-1)*nSteps+1,nframes);
    hiframe=min(nSteps*ichunk,nframes);
    %read the lower frame
    BFFrameRaw = double(read(behavior_vidobj,lowframe));
    BFFrameRaw=BFFrameRaw(:,:,1); %for old videos images come out as RGB, but we only need one
    %select centerline points
    display(['Select Points for frame ' num2str(ichunk) ' of ' num2str(nCells+1) ', showing frame ' num2str(lowframe)]);
    imagesc(BFFrameRaw);
    colormap gray
    tmpS = size(BFFrameRaw);
    hold on
    plot(tmpS(:,1)/2, tmpS(:,2)/2, 'rx')
    set(gcf,'Name',['Select Points for frame ' num2str(ichunk) ])

    [xpts,ypts]=getpts();
    clStartI{ichunk}=[xpts,ypts];
    pause(.3)
    
    if ichunk<=nCells
        bfCell_i{ichunk}=lowframe:hiframe;
    end
    sample_images(:,:,ichunk)=BFFrameRaw;
end
close all
mean_sample=sum(sample_images,3);

%% preview centerlines

close all
button = questdlg('Would you like to preview the centerlines?');
if strcmp(button,'Yes')
    for i=1:nCells+1
        lowframe=min((i-1)*nSteps+1,nframes);
        BFFrameRaw = double(read(behavior_vidobj,lowframe));
        %  BFFrame(BFFrame<0)=0;
        imagesc(BFFrameRaw)
        hold on
        plot(clStartI{i}(:,1),clStartI{i}(:,2),'r');
        pause(.3)
        hold off
    end
end


%% find which pixels of the subset vary the most, excluding ones that may
% bubble related
button = questdlg('Is there an obvious bubble in the video?');
bubble_mask=false(bf_imsize);
while strcmp(button,'Yes')
    display('Crop out where the bubble explores')
    imagesc(mean_sample.*~bubble_mask)
    bubble_mask=and(bubble_mask,roipoly());
    button = questdlg('Is there another bubble?');
end

%% Remove worm for background calculation if it doesnt move much
button = questdlg('Is the worm moving a lot??');
if strcmp(button,'No')
    display('Crop out the worm!')
else
    display('The worms head is always in the center of the image, crop out small head region in the center');
end
imagesc(mean_sample)
worm_mask=roipoly();

%% calculate  background for fluor
progressbar(0);
fluor_stack=0;
skip=max(500,round(nframes_fluor/1000)); %don't need to do every frame, only do 1 every skip.
counter=0;
for itime=1:skip:nframes_fluor;
    progressbar(itime/nframes_fluor);
    fluor_frame = read(fluor_vidobj,itime);
    fluor_stack=fluor_stack+double(fluor_frame(:,:,1));
    counter=counter+1;
    
end
progressbar(1);
fluor_stack_proj=(fluor_stack/counter);
%% pick a region around the head to cut out.
%cut out the head and reinperpolate in to get a background without the
%brain
display('Select an ROI around the bright brain region')
imagesc(fluor_stack_proj);
hold on
fluor_mask=roipoly();
hold off

%% compile masks

masks.fluor_mask=fluor_mask;
masks.worm_mask=worm_mask;
masks.bubble_mask=bubble_mask;

%% save all outputs

display('Done! if you need tips, run wormCL_tip_clicker, otherwise, move the file:')
display([low_mag_folder filesep 'CLworkspace'])
display('into the same folder in /tigress/LEIFER/PanNeuronal')
save([low_mag_folder filesep 'CLworkspace'],...
    'clStartI',...
    'masks',...
    'nCells',...
    'cline_para',...
    'fluormovie',...
    'behaviormovie');

