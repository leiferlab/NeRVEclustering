function clusterCL_start(dataFolder)
%ClusterCL_start takes in an input of a dataFolder and prepares the data
%for submission into the centerline finding program. The dataFolder
%requires the LowMag folder to have a CLworkspace.mat file created by
%initializeCLWorkspace.m. It can also have a tip_coordinates.mat file. 
%It classifies frames into different background states and calculates
%backgrounds for them. It also sets up which frames each job will analyze
%in the subsequent step where clusterWormCenterline is called on multiple
%threads. 



%% write time start stamp
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
if contains(hostname,'della')
    Fid=fopen([dataFolder filesep 'status.txt'],'a');
    status=[datestr(datetime('now')) ': Starting CL initial processing \n'];
    fprintf(Fid,status);
    fclose(Fid);
end




%% find CL workspace with masks and initial centerlines

workspace_file=dir([dataFolder filesep 'CLworkspace*']);

if isempty(workspace_file)
    workspace_file=dir([dataFolder filesep '*' filesep 'CLworkspace*']);
    if isempty(workspace_file)
       error(...
            ['LowMag folder is missing! ensure the Low mag folder is'...
            ' in the BrainScanner Folder and the CLworkspace is present in it'])
    end
    
end

workspace_file=workspace_file(1);
workspace_file=[workspace_file.folder filesep workspace_file.name];
low_mag_folder=fileparts(workspace_file);

%% load the CL workspace
cl_workspace=load(workspace_file);
masks=cl_workspace.masks;
clStartI=cl_workspace.clStartI;
nCells=cl_workspace.nCells;
cline_para=cl_workspace.cline_para;
%% load tip file if present

display('Load tip file if present, otherwise, cancel')
tip_file=[low_mag_folder filesep 'tip_coodinates.mat'];
if exist(tip_file,'file')
    newtips=load(tip_file);
    display(' Tip file found, loading tips');
else
    newtips=[];
    display('No tips found!');
end

%% get video files

%setup paths to movies
fluormovie=cl_workspace.fluormovie;
behaviormovie=cl_workspace.behaviormovie;
%initialize video objects
fluor_vidobj= VideoReader(fluormovie);
behavior_vidobj = VideoReader(behaviormovie);

%get video data
nframes=round(behavior_vidobj.Duration*behavior_vidobj.FrameRate);
nframes_fluors=round(fluor_vidobj.Duration*fluor_vidobj.FrameRate);

bf_imsize=[behavior_vidobj.Height,behavior_vidobj.Width];

%% take subsample of ref pointsfrom each behavior frame
%make ref points
refpoints_x=meshgrid(1:20:bf_imsize(1));
refpoints_y=refpoints_x';
%only take points that dont over lap with the bubbles
mask_points=find(~masks.bubble_mask);
refpoints=sub2ind(bf_imsize,refpoints_x(:),refpoints_y(:));
refpoints=refpoints(ismember(refpoints,mask_points));

%get intensities for all of the reference points in every frame
refintensity=nan(length(refpoints),nframes);
parfor itime=1:nframes
    behavior_vidobj_par = VideoReader(behaviormovie);
    bf_frame =read(behavior_vidobj_par,itime);
    refintensity(:,itime)=bf_frame(refpoints);
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
refintensityZ(refintensityZ>25)=25;
refintensityZ(refintensityZ<-25)=-25;

%
[~,score,~,~] = pca(refintensityZ');

% cluster backgrounds into 11 groups, calculate a background for each of
% them
frame_bg_lvl=kmeans(score(:,1:3),11);

frame_bg_list=unique(frame_bg_lvl);
frame_bg_list=frame_bg_list(~isnan(frame_bg_list));

frames_per_bg=hist(frame_bg_lvl,1:max(frame_bg_lvl));

%if some of them are underpopulated, re-classify them into another cluster
%using simple knn.
while any(frames_per_bg< 200 & frames_per_bg>0)
    target=find(frames_per_bg==min(frames_per_bg));
    trainers=frame_bg_lvl~=target;
    mdl=fitcknn(score(trainers,1:3),frame_bg_lvl(trainers));
    class= predict(mdl,score(~trainers,1:3));
    frame_bg_lvl(~trainers)=class;
    frames_per_bg=hist(frame_bg_lvl,1:max(frame_bg_lvl));

end

%% Now, for each cluster, go through and average all frames from that cluster

mean_bf_all=nan(bf_imsize(1),bf_imsize(2),length(frame_bg_list));

%parfor loop over different background levels and calculate the mean image.
%
parfor i_bg=1:length(frame_bg_list)
    %randomly sample 500 backgrounds, makes things quicker for longer
    %movies
    time_list=find(frame_bg_lvl==frame_bg_list(i_bg));
    time_list=randsample(time_list,min(500,length(time_list)),0);
    mean_bf=zeros(bf_imsize);
    behavior_vidobj_par = VideoReader(behaviormovie);
    
    %add up each frame in a certain level
    for i_time=1:length(time_list)
        currentTime=time_list(i_time);
        bf_frame = read(behavior_vidobj_par,currentTime);
        mean_bf=mean_bf+double(bf_frame(:,:,1));
    end
    %calculate average.
    mean_bf_all(:,:,i_bg)=(mean_bf/length(time_list));
end

% remove worm if worm is present and not moving

for iZ=1:size(mean_bf_all,3)
    temp=mean_bf_all(:,:,iZ);
    temp(masks.worm_mask)=nan;
    temp=inpaint_nans(temp); % interpolate to fill in where worm doesnt move
    mean_bf_all(:,:,iZ)=temp;
end


display('finished!')

%% calculate  background for fluor
%take average
progressbar(0);
fluor_stack=0;
skip=max(20,round(nframes_fluors/1000)); %don't need to do every frame, only do 1 every skip.
counter=0;
for itime=1:skip:nframes_fluors
    progressbar(itime/nframes_fluors);
    fluor_frame = read(fluor_vidobj,itime);
    fluor_stack=fluor_stack+double(fluor_frame(:,:,1));
    counter=counter+1;
    
end
progressbar(1);
f_background=(fluor_stack/counter);

%fill in holes with smoothed interpolation
f_background(masks.fluor_mask)=nan;
f_back_smooth=smooth2a(f_background,15,15);
f_back_smooth=inpaint_nans(f_back_smooth);
f_background(masks.fluor_mask)=mean(f_back_smooth(masks.fluor_mask));

%% make lists of frames each thread will find centerlines for
% flip the bfcell and interleave so that each cell is repeated once
%forward, once back
nSteps=ceil(nframes/nCells); %number of frames in each chunk (except last)
% number of initialized centerlines
bfCell_i=cell(1,nCells);
for ichunk=1:nCells+1
    %get bounds for each chunk
    lowframe=min((ichunk-1)*nSteps+1,nframes);
    hiframe=min(nSteps*ichunk,nframes);
    %make the list go from the lowframe to the high frame
    if ichunk<=nCells
        bfCell_i{ichunk}=lowframe:hiframe;
    end
end

%recreate that cell, but now with the frames counting down
bfCellRev=cellfun(@(x) fliplr(x), bfCell_i, 'uniform',0);

%interleave those two so cell 1 and 2 have the same frames, but with one
%going back and one going forward, etc

bf_list_cell=reshape([bfCell_i;bfCellRev],1,[]);


initial_cl=reshape([clStartI,circshift(clStartI,[0 -1])]',1,[]);
initial_cl=initial_cl(2:end-1);


%% also set reference length based on initializations
worm_length_fun= @(x) sum(sqrt(sum(diff(x).^2,2)));
w_lengths=cellfun(@(x) worm_length_fun(x), initial_cl);
cline_para.refL=mean(w_lengths)/100;


%% Add Fourier mask for low mag behavior image
% Xinwei add for reduce chip posts background.
FourierMask=GenerateFourierMask(behavior_vidobj);
%%

load([low_mag_folder filesep 'CLworkspace']);
%%
% to avoid overwriting tips with empty array if rerunning centerline scipt
tips = newtips;
save([low_mag_folder filesep 'CLworkspace'],...
    'bf_list_cell',...
    'mean_bf_all',...
    'f_background',...
    'initial_cl',...
    'flash_loc_idx',...
    'frame_bg_lvl',...
    'cline_para',...
    'tips',...
    'FourierMask',...
    '-append');




%% write time end stamp
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
if contains(hostname,'della')
    Fid=fopen([dataFolder filesep 'status.txt'],'a');
    status=[datestr(datetime('now')) ': Finished CL initial processing \n'];
    fprintf(Fid,status);
    fclose(Fid);
end

