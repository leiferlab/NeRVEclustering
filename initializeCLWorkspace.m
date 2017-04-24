%% initializeCLWorkspace
% This script is used for initializing parameters and data for centerline
% detection. This script can be run cell by cell, with a few or the cells
% running for or parfor loops over every frame of the video, and some cells
% require user inpurts to initialize centerlines or crop out image regions.
% Prior to running this script, you may run wormCL_tip_clicker.m for a GUI
% that allows the user to click centerline tips. 



%% Initialize fitting parameters for centerline, 

cline_para.refIdx=10;
cline_para.tipRegion=45;
cline_para.endRepulsion=.3;
cline_para.repulsionD=20;
cline_para.heat=3;
cline_para.CLalpha=5;
cline_para.CLbeta=100;
cline_para.gamma=25;  
cline_para.kappa=60;
cline_para.endkappa=5;
cline_para.gradient_force=20;
cline_para.showFlag=0;
cline_para.iterations=400;

cline_para.stretching_force_factor=[.3 .3];
cline_para.refSpring=.01;
cline_para.stretch_ends_flag=1;
cline_para.refL=6;
cline_para.memForce=.05;


nCells=16;

smoothkernal=gausswin(1000, 4);
smoothkernal=smoothkernal/sum(smoothkernal);


%% Select datafolder for analysis
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% get lowmag folder
low_mag_folder=dir([dataFolder filesep 'LowMag*']);
if isempty(low_mag_folder)
    error(...
    'LowMag folder is missing! ensure the Low mag folder is in the BrainScanner Folder')
end
low_mag_folder=[dataFolder filesep low_mag_folder(1).name];


%% load alignment data

alignments=load([dataFolder filesep 'alignments.mat']);
alignments=alignments.alignments;

%% load tip file if present

display('Load tip file if present, otherwise, cancel')
tip_file=[low_mag_folder filesep 'tip_coodinates.mat'];
if exist(tip_file,'file')
    tips=load(tip_file);
    display(' Tip file found, loading tips');
else
    tips=[];
    display('No tips found!');
end

%% set up low magvideos, we've changed the way we save data, the older version
%includes a HUDS avi file and yaml files for metadata

aviFiles=dir([low_mag_folder filesep '*.avi']);
aviFiles={aviFiles.name}';
HUDFiles=aviFiles(cellfun(@(x) ~isempty(strfind(x,'HUDS')),aviFiles));
oldFlag=length(HUDFiles);

if ~oldFlag
    %import textdata
    camdata=importdata([ low_mag_folder  filesep 'CamData.txt']);
    %get timing of each frames
    time=camdata.data(:,2);
    %setup paths to movies
    fluormovie=[low_mag_folder filesep 'cam0.avi'];
    behaviormovie=[low_mag_folder filesep 'cam1.avi'];
    %get movie length
    nframes=length(camdata.data);
    bf2fluor_lookup=[];
else
    aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviormovie=[dataFolder filesep name filesep aviFiles{~aviFluorIdx}];
    fluormovie=[dataFolder filesep name filesep aviFiles{aviFluorIdx}];
     %% set up timing alignments and lookups

    %get timing sync for old data movies, folders were hard saved at
    %1200x600
    [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);
    nframes=length(bfAll.frameTime);
    fluorIdxList=1:length(fluorAll.frameTime);
    bf2fluor_lookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');
    
end



%initialize video objects
behavior_vidobj = VideoReader(behaviormovie);
fluor_vidobj= VideoReader(fluormovie);


%% take subsample of ref pointsfrom each behavior frame
bf_imsize=[behavior_vidobj.Height,behavior_vidobj.Width];
%make ref points
refpoints_x=meshgrid(1:20:bf_imsize(1));
refpoints_y=refpoints_x';
refpoints=sub2ind(bf_imsize,refpoints_x(:),refpoints_y(:));
refintensity=nan(length(refpoints),nframes);
parfor_progress(nframes);
parfor itime=1:nframes;
    % if ~any(fluorAll.flash_loc==itime)
    behavior_vidobj_par = VideoReader(behaviormovie);
    bf_frame =read(behavior_vidobj_par,itime);
    refintensity(:,itime)=bf_frame(refpoints);
    parfor_progress;
    %   end
end

%% find flash locations by looking at intensity of reference points 4sdev above mean
refintensity_mean=mean(refintensity);
refintensityZ_frames=zscore(refintensity_mean);
flash_loc=refintensityZ_frames>4;
flash_loc_idx=find(flash_loc);

%% find which pixels of the subset vary the most, excluding ones that may 
% bubble related
display(['Crop out region with background imformation ( exclude bubbles and obvious worm)'])
refi_std=nanstd(refintensity,[],2);
std_im=reshape(refi_std,55,55);
imagesc(std_im)
std_mask=roipoly();
std_pix=std_im(std_mask);

important_pix= std_mask(:);


%% do PCA on reference point intensities to classify background frames
%get z scores
refintensity_mean=nanmean(refintensity(:,~flash_loc),2);
refintensityZ=bsxfun(@minus, double(refintensity),refintensity_mean);
refintensityZ(:,flash_loc)=0;
refintensityZ=refintensityZ(important_pix,:);
refintensityZ=zscore(refintensityZ,0,2);
%do
[coeff,score,latent,tsquared] = pca(refintensityZ');

% cluster backgrounds into 11 groups, calculate a background for each of
% them
frame_bg_lvl=kmeans(score(:,1:2),11);

frame_bg_list=unique(frame_bg_lvl);
frame_bg_list=frame_bg_list(~isnan(frame_bg_list));
%% calculate multiple background lvls for BF
%progressbar(0);
mean_bf_all=nan(bf_imsize(1),bf_imsize(2),length(frame_bg_list));
parfor_progress(nframes);

%parfor loop over different background levels and calculate the mean image.
%
parfor i_bg=1:length(frame_bg_list);
    time_list=find(frame_bg_lvl==frame_bg_list(i_bg));
    fluor=zeros(bf_imsize);
    behavior_vidobj_par = VideoReader(behaviormovie);
    
    %add up each frame in a certain level
    for i_time=1:length(time_list);
        currentTime=time_list(i_time);
        bf_frame = read(behavior_vidobj_par,currentTime);
        fluor=fluor+double(bf_frame(:,:,1));
        parfor_progress;
    end
    %calculate average.
    mean_bf_all(:,:,i_bg)=(fluor/length(time_list));
end
display('finished!')

%% calculate  background for fluor
progressbar(0);
fluor_stack=0;
skip=round(nframes/1000); %don't need to do every frame, only do 1 every skip.
counter=0;
for itime=1:skip:nframes;
    progressbar(itime/nframes);
    if ~any(flash_loc_idx==itime)
        fluor_frame = read(fluor_vidobj,itime);
        fluor_stack=fluor_stack+double(fluor_frame(:,:,1));
        counter=counter+1;
    end
end
progressbar(1);
fluor_stack_proj=(fluor_stack/counter);
%% pick a region around the head to cut out.
%cut out the head and reinperpolate in to get a background without the
%brain
display('Select an ROI around the bright brain region')

imagesc(fluor_stack_proj);
hold on
head_rectangle=roipoly();
f_background=fluor_stack_proj;
f_background(head_rectangle)=nan;
f_back_smooth=smooth2a(f_background,15,15);
f_back_smooth=inpaint_nans(f_back_smooth);
f_background(head_rectangle)=mean(f_back_smooth(head_rectangle));
imagesc(f_background);
hold off
%% pick a region around the head to cut out, now for behavior
%cut out the head and reinperpolate in to get a background without the
%brain
display('Select an ROI around where the worm is to try and get him out of the background')
imagesc(max(mean_bf_all,[],3));
head_rectangle=roipoly();

for iZ=1:size(mean_bf_all,3);
    
    temp=mean_bf_all(:,:,iZ);
    temp(head_rectangle)=nan;
    temp=inpaint_nans(temp);
    mean_bf_all(:,:,iZ)=temp;
end
meanBfBW=mean_bf_all>90; %flip signs in bright area
imagesc(mean(mean_bf_all,3));
%% make filter background, might not be used yet.
mean_bf_all_filtered=mean_bf_all;
for iZ=1:size(mean_bf_all_filtered,3);
    temp=mean_bf_all_filtered(:,:,iZ);
    temp=bpass(temp,1,40);
    mean_bf_all_filtered(:,:,iZ)=temp;
end

%% initialize centerline points

%cut up the entire video into nCells chunks and then initialize nCells+1
%centerlines for fitting each chunk forwards and backwards.

nSteps=ceil(nframes/nCells); %number of frames in each chunk (except last)
bfCell_i=cell(nCells,1);
clStartI=cell(nCells+1,1);
display(['Click the points of the centerline starting at the head. When '...
    'you are done, double click the last point.']);
%%
for ichunk=1:nCells+1
    %get bounds for each chunk
    lowframe=min((ichunk-1)*nSteps+1,nframes);
    hiframe=min(nSteps*ichunk,nframes);
    %read the lower frame
    BFFrameRaw = double(read(behavior_vidobj,lowframe));

   % background_raw=mean_bf_all(:,:,frame_bg_lvl(lowframe));
    %scale background for best match before subtraction.
   % c=sum(sum(BFFrameRaw.*background_raw))/sum(background_raw(:).^2);
   % BFFrameRaw=BFFrameRaw-background_raw*c;    
    %select centerline points
    display(['Select Points for frame ' num2str(ichunk) ' of ' num2str(nCells+1) ', showing frame ' num2str(lowframe)]);
    imagesc(BFFrameRaw);
    [xpts,ypts]=getpts();
    clStartI{ichunk}=[xpts,ypts];
    pause(.3)
    
    if ichunk<=nCells
        bfCell_i{ichunk}=lowframe:hiframe;
    end
    
end
close all


%% preview centerlines
close all
for i=1:nCells+1
    lowframe=min((i-1)*nSteps+1,nframes);
    BFFrameRaw = double(read(behavior_vidobj,lowframe));
    BFFrame = BFFrameRaw-mean_bf_all(:,:,frame_bg_lvl(lowframe));
  %  BFFrame(BFFrame<0)=0;
    imagesc(BFFrame)
    hold on
    %plot(clStartI{i}(:,1),clStartI{i}(:,2),'r');
    pause(.3)
    hold off
end

%%
% flip the bfcell and interleave so that each cell is repeated once
%forward, once back

bfCellRev=cellfun(@(x) fliplr(x), bfCell_i, 'uniform',0);
bf_list_cell=reshape([bfCell_i,bfCellRev]',1,[]);
initial_cl=reshape([clStartI,circshift(clStartI,[0 -1])]',1,[]);
initial_cl=initial_cl(2:end-1);


%% also set reference length based on initializations
worm_length_fun= @(x) sum(sqrt(sum(diff(x).^2,2)));
w_lengths=cellfun(@(x) worm_length_fun(x), initial_cl);
cline_para.refL=mean(w_lengths)/100;
%%
save([dataFolder filesep 'CLworkspace'],...
    'bf_list_cell',...
    'mean_bf_all',...
    'f_background',...
    'initial_cl',...
    'flash_loc_idx',...
    'frame_bg_lvl',...
    'cline_para',...
    'bf2fluor_lookup',...
    'tips');


