
%% Initialize fitting parameters

maxFrame=Inf; %%% CHECK FOR LAST USEFUL FRAME IN VIDEO

cline_para.refIdx=10;
cline_para.tipRegion=45;
cline_para.endRepulsion=.3;
cline_para.repulsionD=20;
cline_para.heat=3;
cline_para.CLalpha=5;
cline_para.CLbeta=100;
cline_para.gamma=25;
cline_para.kappa=5;
cline_para.endkappa=2;
cline_para.gradient_force=20;
cline_para.showFlag=0;
cline_para.iterations=400;

cline_para.stretching_force_factor=[.3 .3];
cline_para.refSpring=.01;
cline_para.stretch_ends_flag=1;
cline_para.refL=6;
cline_para.memForce=.005;


nCells=16;

smoothkernal=gausswin(1000, 4);
smoothkernal=smoothkernal/sum(smoothkernal);


%% Select datafolder for analysis
dataFolder=uipickfiles();
dataFolder=dataFolder{1};

%% load alignment data
alignmentLocation='Y:\CommunalCode\3dbrain\registration';
backgroundLocation='Y:\CommunalCode\3dbrain\background';
if exist([dataFolder filesep 'alignments.mat'],'file')
    %try to load the alignment file in the folder, otherwise, select them
    %individual in the registration folder
    alignments=load([dataFolder filesep 'alignments.mat']);
    alignments=alignments.alignments;
else
    %Select the alignment files in the order prompted
    display('Select Low Res Alignment')
    
    lowResFluor2BF=uipickfiles('FilterSpec',alignmentLocation);
    lowResFluor2BF=load(lowResFluor2BF{1});
    lowResBF2FluorT=invert(lowResFluor2BF.t_concord);
    
    display('Select Hi to Low Fluor Res Alignment')
    Hi2LowResF=uipickfiles('FilterSpec',alignmentLocation);
    Hi2LowResF=load(Hi2LowResF{1});
    
    display('Select Hi Res Alignment')
    
    S2AHiRes=uipickfiles('FilterSpec',alignmentLocation);
    S2AHiRes=load(S2AHiRes{1});
    rect1=S2AHiRes.rect1;
    rect2=S2AHiRes.rect2;
    
    % if there's a background image, load it as well into alignments.
    display('select a background image for this size himag video');
    
    backgroundImage=uipickfiles('FilterSpec',backgroundLocation);
    if iscell(backgroundImage)
        backgroundImage=load(backgroundImage{1});
        backgroundImage=backgroundImage.backgroundImage;
    else
        backgroundImage=0;
    end
    
    
    %% if you select them individually, bundle them and save it for later use
    alignments.lowResFluor2BF=lowResFluor2BF;
    alignments.S2AHiRes=S2AHiRes;
    alignments.Hi2LowResF=Hi2LowResF;
    alignments.background=backgroundImage;
    save([dataFolder filesep 'alignments'],'alignments');
    
end

%% load tip file if present
try
tip_file=uipickfiles('filterspec', dataFolder);
tip_file=tip_file{1};
tips=load(tip_file);
catch
    tips=[];
end

%% set up low magvideos, we've changed the way we save data, the older version
%includes a HUDS avi file and yaml files for metadata

aviFiles=dir([dataFolder filesep '*.avi']);
aviFiles={aviFiles.name}';
HUDFiles=aviFiles(cellfun(@(x) ~isempty(strfind(x,'HUDS')),aviFiles));
oldFlag=length(HUDFiles);

if ~oldFlag
    
    %import textdata
    camdata=importdata([ dataFolder  filesep 'camData.txt']);
    %get timing of each frames
    time=camdata.data(:,2);
    %setup paths to movies
    fluormovie=[dataFolder filesep 'cam0.avi'];
    behaviormovie=[dataFolder filesep 'cam1.avi'];
    %get movie length
    nframes=length(camdata.data);
    bf2fluor_lookup=[];
else
    aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviormovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluormovie=[dataFolder filesep aviFiles{aviFluorIdx}];
     %% set up timing alignments and lookups

    %get timing sync for old data movies, folders were hard saved at
    %1200x600
    [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[1200 600]);
    nframes=length(bfAll.frameTime);
    fluorIdxList=1:length(fluorAll.frameTime);
    bf2fluor_lookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');
    
end



%initialize video objects
behavior_vidobj = VideoReader(behaviormovie);
fluor_vidobj= VideoReader(fluormovie);


%% calculate backgrounds and projections
bf_imsize=[behavior_vidobj.Height,behavior_vidobj.Width];
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

%% find which pixels of the subset vary the most
refi_std=nanstd(refintensity,[],2);
important_pix=refi_std>quantile(refi_std,.9);


%% do PCA on reference point intensities to classify background frames
%get z scores
refintensity_mean=nanmean(refintensity(:,~flash_loc),2);
refintensityZ=bsxfun(@minus, double(refintensity),refintensity_mean);
refintensityZ(:,flash_loc)=nan;
refintensityZ=refintensityZ(important_pix,:);
%do PCA
[coeff,score,latent,tsquared] = pca(refintensityZ');

%get projection onto first PC for all frames
frame_zscore=colNanFill(score(:,1));

%filter scores
frame_zscore=frame_zscore-conv(frame_zscore,smoothkernal,'same');
%remove flashes
frame_zscore(flash_loc)=nan;
frame_zscore=frame_zscore/nanstd(frame_zscore);

%remove outliers for better quantiles
frame_zscore(abs(frame_zscore)>4)=nan;
%get the quantil each frame falls in,
boundaryI=quantile(frame_zscore,10);
frame_bg_lvl=sum(bsxfun(@le,frame_zscore,boundaryI),2);
frame_bg_lvl=frame_bg_lvl+1;
%testing new way to cluster backgrounds)
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
display('Select an ROI around the bright brain region')
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

    background_raw=mean_bf_all(:,:,frame_bg_lvl(lowframe));
    %scale background for best match before subtraction.
    c=sum(sum(BFFrameRaw.*background_raw))/sum(background_raw(:).^2);
    BFFrameRaw=BFFrameRaw-background_raw*c;    
    %select centerline points
    display(['Select Points for frame ' num2str(ichunk) ' of ' num2str(nCells+1)]);
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
    BFFrameRaw = read(behavior_vidobj,lowframe);
    imagesc(BFFrameRaw)
    hold on
    plot(clStartI{i}(:,1),clStartI{i}(:,2));
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



%%
CLcell=cell(1,2*nCells);
CL_I=CLcell;
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


