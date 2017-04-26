function clusterCL_start(dataFolder)

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

%% get low mag folder
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
%% load CL workspace with masks and centerlines
cl_workspace=load([dataFolder filesep 'CLworkspace']);
masks=cl_workspace.masks;
clStartI=cl_workspace.clStartI;
nCells=cl_workspace.nCells;
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

%%
%import textdata
camdata=importdata([ low_mag_folder  filesep 'CamData.txt']);
%setup paths to movies
fluormovie=[low_mag_folder filesep 'cam0.avi'];
behaviormovie=[low_mag_folder filesep 'cam1.avi'];


%initialize video objects
fluor_vidobj= VideoReader(fluormovie);
behavior_vidobj = VideoReader(behaviormovie);

nframes=round(behavior_vidobj.Duration*behavior_vidobj.FrameRate);
bf_imsize=[behavior_vidobj.Height,behavior_vidobj.Width];

%% take subsample of ref pointsfrom each behavior frame
%make ref points
refpoints_x=meshgrid(1:20:bf_imsize(1));
refpoints_y=refpoints_x';
mask_points=find(~masks.bubble_mask);
refpoints=sub2ind(bf_imsize,refpoints_x(:),refpoints_y(:));
refpoints=refpoints(ismember(refpoints,mask_points));
refintensity=nan(length(refpoints),nframes);
parfor_progress(nframes);
parfor itime=1:nframes;
    behavior_vidobj_par = VideoReader(behaviormovie);
    bf_frame =read(behavior_vidobj_par,itime);
    refintensity(:,itime)=bf_frame(refpoints);
    parfor_progress;
end

%% do PCA on reference point intensities to classify background frames
%get z scores

% find flash locations by looking at intensity of reference points 4sdev above mean
refintensity_mean=mean(refintensity);
refintensityZ_frames=zscore(refintensity_mean);
flash_loc=refintensityZ_frames>4;
flash_loc_idx=find(flash_loc);



refintensity_mean=nanmean(refintensity(:,~flash_loc),2);
refintensityZ=bsxfun(@minus, double(refintensity),refintensity_mean);
refintensityZ(:,flash_loc)=0;
refintensityZ=zscore(refintensityZ,0,2);
%
[coeff,score,latent,tsquared] = pca(refintensityZ');

% cluster backgrounds into 11 groups, calculate a background for each of
% them
frame_bg_lvl=kmeans(score(:,1:2),11);

frame_bg_list=unique(frame_bg_lvl);
frame_bg_list=frame_bg_list(~isnan(frame_bg_list));



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

% remove worm if worm is present

for iZ=1:size(mean_bf_all,3);
    
    temp=mean_bf_all(:,:,iZ);
    temp(masks.worm_mask)=nan;
    temp=inpaint_nans(temp);
    mean_bf_all(:,:,iZ)=temp;
end


display('finished!')

%% for fluor, go through all then take off mask

%% calculate  background for fluor
progressbar(0);
fluor_stack=0;
skip=max(20,round(nframes/1000)); %don't need to do every frame, only do 1 every skip.
counter=0;
for itime=1:skip:nframes;
    progressbar(itime/nframes);
    fluor_frame = read(fluor_vidobj,itime);
    fluor_stack=fluor_stack+double(fluor_frame(:,:,1));
    counter=counter+1;
    
end
progressbar(1);
fluor_stack_proj=(fluor_stack/counter);
f_background=fluor_stack_proj;
f_background(masks.fluor_mask)=nan;
f_back_smooth=smooth2a(f_background,15,15);
f_back_smooth=inpaint_nans(f_back_smooth);
f_background(masks.fluor_mask)=mean(f_back_smooth(masks.fluor_mask));

%%
% flip the bfcell and interleave so that each cell is repeated once
%forward, once back
nSteps=ceil(nframes/nCells); %number of frames in each chunk (except last)
bfCell_i=cell(1,nCells+1);
for ichunk=1:nCells+1
    %get bounds for each chunk
    lowframe=min((ichunk-1)*nSteps+1,nframes);
    hiframe=min(nSteps*ichunk,nframes);
    
    if ichunk<=nCells
        bfCell_i{ichunk}=lowframe:hiframe;
    end
end

bfCellRev=cellfun(@(x) fliplr(x), bfCell_i, 'uniform',0);
bf_list_cell=reshape([bfCell_i,bfCellRev]',1,[]);
initial_cl=reshape([clStartI,circshift(clStartI,[0 -1])]',1,[]);
initial_cl=initial_cl(2:end-1);


%% also set reference length based on initializations
worm_length_fun= @(x) sum(sqrt(sum(diff(x).^2,2)));
w_lengths=cellfun(@(x) worm_length_fun(x), initial_cl);
cline_para.refL=mean(w_lengths)/100;
%%
save([low_mag_folder filesep 'CLworkspace'],...
    'bf_list_cell',...
    'mean_bf_all',...
    'f_background',...
    'initial_cl',...
    'flash_loc_idx',...
    'frame_bg_lvl',...
    'cline_para',...
    'tips',...
    '-append');


