function clusterWormCenterline(dataFolder,iCell,isChip,show2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clusterworm tracker fits centerlines to behavior videos from our whole
%brain imaging setup. a CL workspace must be loaded with initial parameters
%and paths in order to run this code, and the activeContourFit program
%requires the eigenworms to be loaded as "eigbasis" on to the main window.
%It also requires an alignment.mat file in the BrainScanner folder.
if nargin<3
    isChip=0;
end

if nargin<4
    show2=0;
end

%% initialize some kernels and flags
STARTFLAG=1;

GAUSSFILTER=fspecial('gaussian',30,5);
GAUSSFILTER2=fspecial('gaussian',50,15);

cm_fluor=[26 26];
sdev_nhood=getnhood(strel('disk',5));

se=strel('disk',11);
se=se.Neighborhood;
k=fspecial('Gaussian',5,2.5);


%% load eigenworms and set as global
eigbasis=load('eigenWorms_full.mat');
eigbasis=eigbasis.eigvecs;
setappdata(0,'eigbasis',eigbasis);

%% get folders that have the avi files

%% find CL workspace with masks and initial centerlines

workspace_file=dir([dataFolder filesep 'CLworkspace*']);

if isempty(workspace_file)
    workspace_file=dir([dataFolder filesep '*' filesep 'CLworkspace*']);
    if isempty(workspace_file)
        error(...
            'LowMag folder is missing! ensure the Low mag folder is in the BrainScanner Folder')
    end
end

workspace_file=workspace_file(1);
workspace_file=[workspace_file.folder filesep workspace_file.name];
low_mag_folder=fileparts(workspace_file);
%%


%% load initial variables from CLworkspace

try
    CLworkspace=load(workspace_file);
    bf_list_cell=CLworkspace.bf_list_cell; %bf_list_cell bfCell
    mean_bf_all=CLworkspace.mean_bf_all; %mean_bf_all meanBfAll2
    f_background=CLworkspace.f_background;  %f_background fluorBackground
    initial_cl=CLworkspace.initial_cl; % initial_cl clStart
    flash_loc_idx=CLworkspace.flash_loc_idx; %flash_loc flashLoc
    frame_bg_lvl=CLworkspace.frame_bg_lvl; %frame_bg_lvl newZ2
    if isChip
        FourierMask=CLworkspace.FourierMask;
    end
    %setup paths to movies
    fluormovie=CLworkspace.fluormovie;
    behaviormovie=CLworkspace.behaviormovie;
    
    
    %set up centerline fitting parameters
    cline_para=CLworkspace.cline_para;
catch me
    error('Files are missing from the CLworkspace, did you run clusterCL_start?')
end


if isfield(CLworkspace,'bf2fluor_lookup')
    bf2fluor_lookup=CLworkspace.bf2fluor_lookup; %for old setup of avi files
else
    bf2fluor_lookup=[]; %otherwise, behavior and fluor videos are in sync
end

refIdx=cline_para.refIdx;
cline_para.showFlag=00;

framelist=bf_list_cell{iCell};
cl_all=zeros(100,2,length(framelist)); %cl_all CLall
cl_intensities=zeros(100,length(framelist)); %cl_intensities IsAll


%load tip data and interpolate missing values
if isfield(CLworkspace,'tips')
    if isfield(CLworkspace.tips,'head_pts')
        head_pts=CLworkspace.tips.head_pts;
        tail_pts=CLworkspace.tips.tail_pts;
        
        tip_mat_length=min(size(head_pts,1),size(tail_pts,1));
        head_pts=head_pts(1:tip_mat_length,:);
        tail_pts=tail_pts(1:tip_mat_length,:);
        
        same_pt=all(head_pts==tail_pts,2);
        
        head_time=find(head_pts(:,1) & ~same_pt);
        head_pts_sub=head_pts(head_time,:);
        
        tail_time=find(tail_pts(:,1) & ~same_pt);
        tail_pts_sub=tail_pts(tail_time,:);
        
        if isempty(head_pts_sub)
            head_pt_list=[];
            display('WARNING: tip file present, but no head points are Found!')
        else
            head_pt_list=interp1(head_time,head_pts_sub,framelist,'pchip');
        end
        
        if isempty(tail_pts_sub)
            tail_pt_list=[];
            display('WARNING: tip file present, but no tail points are Found!')
        else
            tail_pt_list=interp1(tail_time,tail_pts_sub,framelist,'pchip');
        end
        
        
        cline_para.stretch_ends_flag=0;
        cline_para.stretching_force_factor=[0 0];
        cline_para.endkappa=0;
        cline_para.endRepulsion=.003;
        cline_para.refSpring=0.01; %turn down pulling to center if tips are clicked
    else
        head_pt_list=[];
        tail_pt_list=[];
    end
else
    head_pt_list=[];
    tail_pt_list=[];
    
end



%% select files video files and load avi files
behavior_vidobj = VideoReader(behaviormovie);
fluor_vidobj= VideoReader(fluormovie);

%% create output folder
outputFolder=[low_mag_folder filesep 'CL_files'];
if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end


%% load alignments
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
lowResFluor2BF=alignments.lowResFluor2BF;

%% main loop
for iframe=1:length(framelist)
    %%
    %get frame from framelist
    itime=framelist(iframe);
    if any(head_pt_list)
        cline_para.head_pt=head_pt_list(iframe,:);
        cline_para.tail_pt=tail_pt_list(iframe,:);
    else
        cline_para.head_pt=[];
        cline_para.tail_pt=[];
    end
    
    tic
    try
        if ~isnan(frame_bg_lvl(itime)) && ~any(flash_loc_idx==itime) && itime>=1
            
            if isChip
                bf_frame_raw = read(behavior_vidobj,itime);
                bf_frame_raw=double(bf_frame_raw(:,:,1));
                % fourier mask filter(if not empty)
                %low pass filter to deal with background.
                lowpassMask=false(size(bf_frame_raw,1),size(bf_frame_raw,2));
                lowpassMask(floor(size(bf_frame_raw,1)/2),floor(size(bf_frame_raw,2)/2))=true;
                se_lowpass=strel('square',5);
                lowpassMask=imdilate(lowpassMask,se_lowpass);
                if ~isempty(FourierMask)
                    FourierMask=FourierMask | lowpassMask;
                else
                    FourierMask=lowpassMask;
                end  
                
                fftimg=fftshift(fft2(bf_frame_raw));
                fftimg(logical(FourierMask))=0;
                bf_frame_raw=real( ifft2( ifftshift(fftimg) ) );
                bf_frame=imfilter(bf_frame_raw,k);
            
                % afew filter steps
                bf_frame_std=stdfilt(bf_frame_raw,sdev_nhood);
                bf_frame_std=normalizeRange(bf_frame_std);
            else
                %% filter behavior images
                bf_frame_raw = read(behavior_vidobj,itime);
                bf_frame_raw=double(bf_frame_raw(:,:,1));
                background_raw=mean_bf_all(:,:,frame_bg_lvl(itime));
                %scale background for best match before subtraction.
                c=sum(sum(bf_frame_raw.*background_raw))/sum(background_raw(:).^2);
                bf_frame_raw=bf_frame_raw-background_raw*c;
                
                
                
                
                bg=normalizeRange(background_raw);
                bg_mask=bg>.5;
                bg_mask=AreaFilter(bg_mask,5000,[],8);
                bg_mask=imclose(bg_mask,true(12));
                bg_mask=imdilate(bg_mask,true(25));
                
                
                C2=imfilter(bf_frame_raw,k);
                [H,D]=hessianMatrix(C2,3);
                [Heig,HeigVec]=hessianEig(H);
                
                
                bf_frame_raw(bf_frame_raw<0)=0;
                
                
                H1=stdfilt(HeigVec{1,1}.*Heig(:,:,1),se);
                H2=stdfilt(HeigVec{2,1}.*Heig(:,:,1),se);
                H=sqrt(H1.^2+H2.^2);
                
                H=H*1000;
                H=pedistalSubtract(H);
                H=H-2*median(H(:));
                H(H<0)=0;
                bf_frame=bf_frame_raw*.7+H*.3;
                
                
                % afew filter steps
                bf_frame_std=stdfilt(bf_frame_raw,sdev_nhood);
                bf_frame_std=normalizeRange(bf_frame_std);
            end
            
            %% filter fluor images
            % this is empty for new setup, not empty for old
            if ~isempty(bf2fluor_lookup)
                fluor_time=round(bf2fluor_lookup(itime));
                %band aid solution for if time's dont align well (behavior
                %starts after fluor), this doesnt happen for new setup
                if isnan(fluor_time) || fluor_time<1;
                    fluor_time=1;
                end
            else
                fluor_time= itime;
            end
            
            %% load fluorescent images
            fluor_frame_raw=read(fluor_vidobj,fluor_time);
            fluor_frame_raw=double(fluor_frame_raw(:,:,1))-f_background;
            fluor_frame_raw(fluor_frame_raw<0)=0;
            fluor_frame=fluor_frame_raw;
            fluor_frame=bpass(fluor_frame,1,40);
            
            %% calculate centroid of fluor image
            fluor_mask=false(size(fluor_frame));
            fluor_mask(round(cm_fluor(2))+(-25:25),...
                round(cm_fluor(1))+(-25:25))=true;
            fluorBW=(fluor_frame>(max(fluor_frame(:))/5));
            fluor_mask=fluorBW.*fluor_mask;
            if (any(fluor_mask(:))) && ~STARTFLAG
                fluorBW=fluor_mask;
            end
            fluorBW=AreaFilter(fluorBW); %select only largest obj
            [cmy,cmx]=find(fluorBW);
            cm_fluor=mean([cmx cmy]);
            cm=transformPointsForward(lowResFluor2BF.t_concord,cm_fluor);
            %% make input image
            
            %gassian filter
            if isChip
                inputImage=bf_frame;
            else
                inputImage=filter2(GAUSSFILTER,bf_frame,'same');
            end
            
            bf_frame_std=normalizeRange(imfilter(bf_frame_std,GAUSSFILTER2));
            bfFrameMask=(bf_frame_std>min(graythresh(bf_frame_std),.7));
            
            %make tip image, to help provide forces of extenstion for tip
            tip_image=bfFrameMask.*inputImage;
            bfFrameLTmask=bfFrameMask & tip_image<bf_frame_std;
            tip_image(bfFrameLTmask)=bf_frame_std(bfFrameLTmask); 
            %% fit centerlines using previous centerline as starting point
            %ActiveContourFit_wormRef4 does all the work.
            if STARTFLAG
                % for first round, just load up initialized CL points as
                % past frame
                cl_old=initial_cl{iCell};
                cl_old=distanceInterp(cl_old,100);
                [cl,Is]=ActiveContourFit_wormRef4(...
                    inputImage,tip_image, cline_para, cl_old,refIdx,cm);
                STARTFLAG=0;
            else
                oldTime=iframe-1;
                while isnan(frame_bg_lvl(framelist(oldTime)));
                    oldTime=oldTime-1;
                end
                cl_old=cl_all(:,:,oldTime);
                [cl,Is,Eout]=ActiveContourFit_wormRef4(...
                    inputImage,tip_image, cline_para, cl_old,refIdx,cm);
                
                cl_dist=mean(sqrt(sum((cl-cl_old).^2,2)));
                %don't let the centerlines go too far
                if cl_dist>25 && isempty(cline_para.head_pt)
                    cl=cl_old;
                end
            end
            
            %% plot some of the results if show2 is 1
            if ~mod(itime,show2)
                subplot(1,2,1)
                imagesc(bf_frame_raw);
                hold on
                plot(cl(:,1),cl(:,2),'r');
                plot(cl([1 end],1),cl([1 end],2),'og');
                scatter(cm(1),cm(2),'gx');
                plot([cl(refIdx,1) cm(1)],[cl(refIdx,2) cm(2)],'g');
                hold off
                subplot(1,2,2);
                imagesc(bf_frame);
                hold on
                plot(cl(:,1),cl(:,2),'r');
                plot(cl([1 end],1),cl([1 end],2),'og');
                scatter(cm(1),cm(2),'gx');
                plot([cl(refIdx,1) cm(1)],[cl(refIdx,2) cm(2)],'g');
                hold off
                drawnow
            end
            
            %% reinterpolate to 100 points
            cl=distanceInterp(cl,100);
            cl_all(:,:,iframe)=cl;
            cl_intensities(:,iframe)=Is;
        else
            %if theres a frame error, reuse the previous frame's CL.
            cl_all(:,:,iframe)=cl_all(:,:,iframe-1);
            cl_intensities(:,iframe)=cl_intensities(:,iframe-1);
        end
        display(['Completed frame '  num2str(itime) ', cell '...
            num2str(iCell) ' in ' num2str(toc) ' s'])
    catch me
        %if error, print error, reuse previous centerline
        display(['error frame ' num2str(itime) ', cell ' num2str(iCell)])
        if cline_para.showFlag==0 && show2==0
            me
        end
        if iframe>1
            cl_all(:,:,iframe)=cl_all(:,:,iframe-1);
            cl_intensities(:,iframe)=cl_intensities(:,iframe-1);
        else
            %if this is the first frame of the list, zero the outputs
            cl_old=initial_cl{iCell};
            cl_old=distanceInterp(cl_old,100);
            cl_all(:,:,iframe)=cl_old;
            cl_intensities(:,iframe)=0;
        end
        
        
        
    end
end
%% save output
outputFilename=[outputFolder filesep 'CL_' num2str(iCell,'%3.2d')];
save(outputFilename,'cl_all','cl_intensities','framelist');
%note: cellList changed to frame list
