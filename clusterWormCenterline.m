function clusterWormCenterline(dataFolder,iCell,show2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clusterworm tracker fits centerlines to behavior videos from our whole
%brain imaging setup. a CL workspace must be loaded with initial parameters
%and paths in order to run this code, and the activeContourFit program
%requires the eigenworms to be loaded as "eigbasis" on to the main window.
%It also requires an alignment.mat file in the BrainScanner folder.

if nargin==2
    show2=0;
end

%% initialize some kernels and flags
STARTFLAG=1;

GAUSSFILTER=fspecial('gaussian',30,5);
GAUSSFILTER2=fspecial('gaussian',50,15);

cm_fluor=[26 26];
sdev_nhood=getnhood(strel('disk',5));


%% load eigenworms and set as global
eigbasis=load('eigenWorms_full.mat');
eigbasis=eigbasis.eigvecs;
setappdata(0,'eigbasis',eigbasis);

%% get folders that have the avi files
%Lowmag folder only exist for new set
d= dir([dataFolder filesep 'LowMagBrain*']);
if isempty(d)
    %for old
    aviFolder=dataFolder;
else
    %for new
    aviFolder=[dataFolder filesep d(1).name];
end
display(aviFolder);

%% load initial variables from CLworkspace
workSpaceFile=[aviFolder filesep 'CLworkspace.mat'];



CLworkspace=load(workSpaceFile);
bf_list_cell=CLworkspace.bf_list_cell; %bf_list_cell bfCell
mean_bf_all=CLworkspace.mean_bf_all; %mean_bf_all meanBfAll2
f_background=CLworkspace.f_background;  %f_background fluorBackground
initial_cl=CLworkspace.initial_cl; % initial_cl clStart
flash_loc_idx=CLworkspace.flash_loc_idx; %flash_loc flashLoc
frame_bg_lvl=CLworkspace.frame_bg_lvl; %frame_bg_lvl newZ2
cline_para=CLworkspace.cline_para;
if isfield(CLworkspace,'bf2fluor_lookup')
    bf2fluor_lookup=CLworkspace.bf2fluor_lookup; %for old setup of avi files
else
    bf2fluor_lookup=[];
end
refIdx=cline_para.refIdx;
cline_para.showFlag=00;

framelist=bf_list_cell{iCell};
cl_all=zeros(100,2,length(framelist)); %cl_all CLall
cl_intensities=zeros(100,length(framelist)); %cl_intensities IsAll


if isfield(CLworkspace,'tips')
    if isfield(CLworkspace.tips,'head_pts')
    head_pts=CLworkspace.tips.head_pts;
    tail_pts=CLworkspace.tips.tail_pts;
    same_pt=all(head_pts==tail_pts,2);
    

    head_time=find(head_pts(:,1) & ~same_pt);
    head_pts_sub=head_pts(head_time,:);
    head_pt_list=interp1(head_time,head_pts_sub,framelist,'pchip');
    
    tail_time=find(tail_pts(:,1) & ~same_pt);
    tail_pts_sub=tail_pts(tail_time,:); 
    tail_pt_list=interp1(tail_time,tail_pts_sub,framelist,'pchip');
    
    cline_para.stretch_ends_flag=0;
    cline_para.endkappa=0;
    
    else
    head_pt_list=[];
    tail_pt_list=[];  
    end
else
    head_pt_list=[];
    tail_pt_list=[];
    
end



%% select files video files and load avi files
if isempty(bf2fluor_lookup)
    %for the newer setup, avi files are in seperate folder
    fluor_movie_file=[aviFolder filesep 'cam0.avi'];
    behavior_movie_file=[aviFolder filesep 'cam1.avi'];
else
    %for older setup
    aviFiles=dir([aviFolder filesep '*.avi']);
    aviFiles={aviFiles.name}';
    aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behavior_movie_file=[aviFolder filesep aviFiles{~aviFluorIdx}];
    fluor_movie_file=[aviFolder filesep aviFiles{aviFluorIdx}];
end

behavior_vidobj = VideoReader(behavior_movie_file);
fluor_vidobj= VideoReader(fluor_movie_file);

%% create output folder
outputFolder=[aviFolder filesep 'CL_files'];
if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end


%% load alignments
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
lowResFluor2BF=alignments.lowResFluor2BF;

%% main loop
for iframe=1:length(framelist);
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
            %% filter behavior images
            bf_frame_raw = read(behavior_vidobj,itime);
            bf_frame_raw=double(bf_frame_raw(:,:,1));
            background_raw=mean_bf_all(:,:,frame_bg_lvl(itime));
            %scale background for best match before subtraction.
            c=sum(sum(bf_frame_raw.*background_raw))/sum(background_raw(:).^2);
            bf_frame_raw=bf_frame_raw-background_raw*c;
            bf_frame_raw=abs(bf_frame_raw);
            % afew filter steps
            bf_frame=imtophat(bf_frame_raw,strel('disk',50));
            bf_frame_std=stdfilt(bf_frame_raw,sdev_nhood);
            bf_frame_std=normalizeRange(bf_frame_std);
            bfstdthresh=(bf_frame_std>graythresh(bf_frame_std));
            bf_frame(bfstdthresh)=abs(bf_frame(bfstdthresh));
            bf_frame_std=bpass((bf_frame_std),1,80);
            bf_frame=bpass(bf_frame,4,80)/10;
            
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
            %combine filtered image and filtered sdev image
            inputImage=2*bf_frame/max(bf_frame(:))+normalizeRange(bf_frame_std);
            inputImage(inputImage<0)=0;
            %gassian filter
            inputImage=filter2(GAUSSFILTER,inputImage,'same');
            inputImage=normalizeRange(inputImage);
            
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
