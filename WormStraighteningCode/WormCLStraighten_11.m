
function [V,pointStats,Vproj,side,xyOffset2,wormBW2]=WormCLStraighten_11...
    (dataFolder,destination,vidInfo,alignments,Vtemplate,zOffset,iStack,side,show)

% takes data for whole brain imaging set, centerlines, himag movies, lowmag
% behavior videos to create straighten worm in given frame
%adding second part to straighten.

%Ver 11, trying to make straightening better with some feedback

%% initial parameters
%size to search around centerline
outputRadius=83.5; %xy search radius around centerline
outputRadiusZ=63.5; %z serach radius around centerilne
outputLength=200; % serach radius along centerline about center

%buffer to add around search radius, to be removed  after alignments
outputRadiusBuff=30;    
outputLengthBuff=100;


%ratio between xypixels and z slices
zRatio=1/3;

%options for segmentation
options.method='invdist';
options.radius=20;
options.power=1;
options.thresh1=.05;
options.minObjSize=50;
options.maxObjSize=400;
options.minSphericity=.80;
options.filterSize=[10 10 4];
options.power=1;
options.prefilter=1;
options.hthresh=0;

%set up destination folder
destination_path=[dataFolder filesep destination];
%% set up different kernals
%gaussians
gaussKernal2=gausswin(200);
gaussKernal2=convnfft(gaussKernal2,gaussKernal2');

%mexican hat filters
Sfilter=max(gaussKernal2(:))-gaussKernal2;
Sfilter(Sfilter<.1)=-(.1-Sfilter(Sfilter<.1))*80;

Sfilter(Sfilter>.8)=0;
Sfilter(Sfilter>0)=1;

Sfilter2=max(gaussKernal2(:))-gaussKernal2;
Sfilter2(Sfilter2<.01)=Sfilter2(Sfilter2<.01)-.3;
Sfilter2(Sfilter2>.6)=0;

Sfilter2(Sfilter2>0)=nnz(Sfilter2<0)/nnz(Sfilter2>0);
Sfilter2=Sfilter2-mean(Sfilter2(:));


%% recover alignments
lowResFluor2BF=alignments.lowResFluor2BF;
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
rect1=S2AHiRes.rect1;

%% set up low fluor video to assist centerline centering

aviFiles=dir([dataFolder filesep '20*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));

d= dir([dataFolder filesep 'LowMagBrain*']);
if ~isempty(d)
    aviFolder=[dataFolder filesep d(1).name];
end
if length(aviFiles)==2
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    fluorMovie=[dataFolder filesep aviFiles{aviFluorIdx}];
elseif isdir(aviFolder)
    fluorMovie=[aviFolder filesep 'cam0.avi'];
else
    
    display('Select avi files,low mag fluor');
    movies=uipickfiles('FilterSpec',dataFolder);
    fluorMovie=movies{1};
end

fluorVidObj= VideoReader(fluorMovie);

%% set up high mag videos
if isempty(vidInfo)
    
    [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);
else
    bfAll=vidInfo.bfAll;
    fluorAll=vidInfo.fluorAll;
    hiResData=vidInfo.hiResData;
    
end
%% set up timing alignments and lookups using interpolation
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(...
    bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(...
    fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');


%% load centerline
%get behavior folder
behaviorFolder=dir([dataFolder filesep 'Behavior*']);
behaviorFolder=behaviorFolder([behaviorFolder.isdir]);
behaviorFolder=[dataFolder filesep behaviorFolder(1).name];

%in behavior folder, get centerline file
centerlineFile=dir([behaviorFolder filesep 'center*']);
centerlineFile=[behaviorFolder filesep centerlineFile(1).name];

%load centerline file and make variables based on field names
centerline=load(centerlineFile);
CLfieldNames=fieldnames(centerline);
CLfieldIdx=cellfun(@(x) ~isempty(strfind(x,'centerline')),CLfieldNames);
CLoffsetIdx=cellfun(@(x) ~isempty(strfind(x,'off')),CLfieldNames);
if any(CLoffsetIdx)
    CLoffset=centerline.(CLfieldNames{CLoffsetIdx});
else
    CLoffset=0;
end
centerline=centerline.(CLfieldNames{CLfieldIdx});


%%
%try the entire pipeline, workspace in case of error for further
%troubleshooting. 
try
    %% load images
    
    datFileDir=dir([dataFolder filesep 'sCMOS_Frames_U16_*']);
    datFile=[dataFolder filesep datFileDir.name];
    [rows,cols]=getdatdimensions(datFile);
    nPix=rows*cols;
    Fid=fopen(datFile);
    
    %select frames to analyze
    hiResIdx=find(hiResData.stackIdx==iStack)+ zOffset;
    %get z values of those frames
    zRange=hiResData.Z(hiResIdx-zOffset);
    zSize=length(hiResIdx);

    %get correspoinding fluor indices
    fluorIdx=round(fluorIdxLookup(hiResIdx));
    fluorIdxRange=[min(fluorIdx) max(fluorIdx)];
    
    %load up lowmag fluor image
    fluor_raw=read(fluorVidObj,fluorIdxRange);
    fluor_raw=squeeze(fluor_raw(:,:,1,:));
    %warp it to align with high mag, it warps to the uncropped high mag so
    %we need to crop it
    fluor_trans=imwarp(fluor_raw,Hi2LowResF.t_concord,...
        'OutputView',Hi2LowResF.Rsegment);
    fluor_trans=fluor_trans((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);

        %do something with status errors!
    status=fseek(Fid,2*(hiResIdx(1))*nPix,-1);
    
    pixelValues=fread(Fid,nPix*(length(hiResIdx)),'uint16',0,'l');
    hiResImage=reshape(pixelValues,rows,cols,length(hiResIdx));
    %subtract background
    hiResImage=bsxfun(@minus, hiResImage,alignments.background);
    
    %% crop and align hi mag images
    worm_im=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    raw_size=size(worm_im);
    %% filter to find center Z (some old things here)
    %3d bpass filter
    worm_smooth=bpass3(worm_im,.5,[20 20 3]);

    %% find middle plane by looking for large thresholded objects. 
    % additional boxcar smoothing in XY
    worm_smooth_thresh=smooth3(worm_smooth,'box',[15,15,1]);
    % thresholded smooth image
    worm_smooth_thresh=normalizeRange(worm_smooth_thresh)>.05;
    
    max_dist_plot=zeros(1,zSize);
    for i=1:zSize
        %for each slice, find object boundary pixels
        [midZx,midZy]=find(bwmorph(worm_smooth_thresh(:,:,i),'remove'));
        %find max distance between those points, the slice with the largest
        %chord is the middle
        max_thresh_dist=max(pdist([midZx midZy]));
        if ~isempty(max_thresh_dist)
            max_dist_plot(i)=max_thresh_dist;
        end
    end
    %find peak in plot, corresponding with slice with largest object
    max_dist_plot=normalizeRange(max_dist_plot(:));
    [~,midZ]=max(smooth(max_dist_plot,3));
    %incase there are repeat vallues, just take the mean (rare)
    midZ=(mean(midZ));
    %round the middle plane in the direction of the middle 
    if midZ>zSize/2
        midZ=floor(midZ);
    else
        midZ=ceil(midZ);
    end
    %dont let middle plane be one of the end planes.
    if midZ==1
        midZ=2;
    elseif midZ>=zSize;
        midZ=zSize-1;
    end

    %% try fix centerline alignemnt by looking at lowmag fluor and finding 
    % a correction offset between the transformed lowmag fluor and the high
    % mag segmentation. 
    
    % get lowmag fluor and threshold, and project
    % cast fluor image
    fluor_trans=normalizeRange(double(fluor_trans));
    %apply automatic threshold and z project
    fluor_thresh=(fluor_trans>(graythresh(fluor_trans(:))));
    fluor_thresh_proj=normalizeRange(sum(fluor_thresh,3));
    
    %remove boundary 
    boundaryRegion=ones(size(fluor_thresh_proj)); 
    boundaryRegion(51:end-50,51:end-50)=0;
    fluor_thresh_proj=fluor_thresh_proj.*~boundaryRegion;
    
 
    %% get proper centerlines to use that correspond to the hiresframes
    hires_range=min(hiResIdx):max(hiResIdx);
    CLIdx=round(bfIdxLookup(hires_range));
    CLIdx=CLIdx-CLoffset;
    % find unique terms, ia are the unique idx, ic is the mapping from the
    % output to the intput.
    [CLIdx,ia,ic]=unique(CLIdx);
    CL_lo=centerline(:,:,CLIdx); %CL in lowmag coordinate
    %% try to fix obvious centerline errors if they're off
    %compare centerlines to median filtered centerelines, replace them if
    %they are far off
    % reshape into 2day array for median filtering
    clx=squeeze(CL_lo(:,1,:));
    cly=squeeze(CL_lo(:,2,:));
    clx_med=medfilt2(clx,[1,5]);
    cly_med=medfilt2(cly,[1,5]);
    
    % cl err is the squared distance from original CL and med filtered 
    cl_err=sqrt(mean((clx_med-clx).^2+(cly_med-cly).^2));
    cl_replace_idx=find(cl_err'>30 & ia>1 & ia<max(ia));
    % if distance is large, replace CL with medfiltered one. 
    CL_lo(:,1,cl_replace_idx)=clx_med(:,cl_replace_idx);
    CL_lo(:,2,cl_replace_idx)=cly_med(:,cl_replace_idx);
    
    
    %% transform centerlines into himag coordinates
    
    %get the centerlines, one for each hires frame, with some repeats
    
    %CL_hi is CL hi mag coordinate, with low mag spacing
    CL_hi=CL_lo(:,:,ic); 
    
    %centerline from behavior coordinate system to fluor coordinate system
    [CL_fluor(:,2,:),CL_fluor(:,1,:)]=transformPointsInverse(...
        lowResFluor2BF.t_concord,CL_lo(:,2,:),CL_lo(:,1,:));
    %centerline from lowmag fluor coordinate system to himag
    [CL_hi(:,1,:),CL_hi(:,2,:)]=transformPointsForward(....
        Hi2LowResF.t_concord,CL_fluor(:,2,:),CL_fluor(:,1,:));
    %subtract the rectangle size
        CL_hi(:,2,:)=CL_hi(:,2,:)-(rect1(2)-1);
    CL_size=size(CL_hi);
    
    %% interpolate to parameterize by length
    %range to interpolate from tip, in highres pixels. This can be generous
    %to start, as the output length is determined byt outputLength. 
    CL_range=2500;
    CL2_hi=zeros(CL_range,2,CL_size(3));
    
    %reinterpolate centerline by length, with spacing appropriate for
    %higher mag coordinate system
    for iCL=1:CL_size(3)
        CL2temp=CL_hi(:,:,iCL);
        %distance steps
        ds=squeeze(sqrt(sum((diff(CL2temp,[],1)).^2,2)));
        s=[0;cumsum(ds) ];
        %reinterpolate from initial distances to pixel distances
        CL2temp=interp1(s,CL2temp,1:1:CL_range,'linear','extrap');
        CL2temp(:,:,iCL)=CL2temp;
        if show
            if iCL==1
                close all
                imagesc(worm_im(:,:,midZ))
                hold on
            end
            plot(CL2temp(:,1),CL2temp(:,2))
            scatter(CL2temp(1500,1),CL2temp(1500,2));
            drawnow
        end
        
    end
    
    %% align centerlines parameterizations by correlation
    % centerlines points can be offset to each other because nothing stops
    % cls from sliding if the tips are not fixed. We need to align the
    % consecutive centerlines to each other.  
    % Alignment of  centerlines using plotMatch2 (distance based) to account
    % for centerline sliding.
    
    %how much to shift each centerline
    shiftVec=zeros(1,length(ia)-1);
    for i=2:length(ia);
        %find best offset between consecutive centerlines.
        CL1temp=CL2_hi(:,:,ia(i));
        CL2temp=CL2_hi(:,:,ia(i-1));
        [corrtemp,r]=plotMatch2(CL1temp,CL2temp);
        [~,shift]=min(corrtemp);
        shiftVec(i)=r(shift);
    end
    % add the offsets back
    shiftVec=cumsum(shiftVec);
    shiftVec=shiftVec-shiftVec(round(length(shiftVec)/2));
    
    %add a buffer in the size of the centerline to allow for shifting
    shift_buffer=500;
    
    %CL2_hi_long is CL2_hi with extra buffer space, buffer cut off later
    CL2_hi_long=zeros(CL_range+shift_buffer+1,2,CL_size(3));
    
    %apply shifts by interpolating
    for iCL=1:CL_size(3)
        CL_temp=CL2_hi(:,:,iCL);
        %also remove shift_buffer.
        new_cl_range=shiftVec(ic(iCL))+(-shift_buffer:CL_range);
        CL_temp=interp1(CL_temp,new_cl_range,'linear','extrap');
        CL2_hi_long(:,:,iCL)=CL_temp;
        if show
            if iCL==1
                close all
                imagesc(worm_im(:,:,midZ))
                hold on
            end
            plot(CL_temp(:,1),CL_temp(:,2))
            scatter(CL_temp(1300,1),CL_temp(1300,2),'w');
            drawnow
        end
    end
    
    %% find offset to better align centerline with hi mag worm brain.
 
    %this occurs in three steps, we find the offset between the brightest
    %centerline point and the brightest image point in lowmag fluor
    %the we find the offset between the centerlines
    %and the lowmag fluorescent images (after transforming into hi mag
    %coordinates), and then finding the offset between the transformed
    %lowmag fluorecent images and the hi mag images. The transformation
    %from bead selection is often off by a bit and it changes over time. If
    %I find a better straightening algorithm then I can dump all of this. 

    [fcm_x, fcm_y]=find(fluor_thresh_proj==max(fluor_thresh_proj));
    fcm_x=mean(fcm_x);
    fcm_y=mean(fcm_y);
    
    %find CL points that has brightest points in low mag fluor
    CL2f_mean=(CL_fluor(:,:,1));
    fluor_raw_proj=sum(fluor_raw,3);
    CL2f_I=interp2(fluor_raw_proj,CL2f_mean(:,2),CL2f_mean(:,1));
    [~,refIdx]=max(CL2f_I);
    %create initial offset as the shift between the high mag fluor center and
    %the brightest index and the high mag CL point found using lowmag data
    xyOffset2=[fcm_y fcm_x]-[CL_hi(refIdx,1,1) CL_hi(refIdx,2,1)];
    
    CL2X=reshape(mean(CL2_hi_long(:,1,1:end),3),[],1,1);
    CL2Y=reshape(mean(CL2_hi_long(:,2,1:end),3),[],1,1);

    %shift around the centerline to optimize overlap between centerline and
    %fluorescence image
    
    %apply filter that emphasizes center of worm brain
    fluor_filt=convnfft(fluor_thresh_proj,Sfilter2,'same');
    % ignore boundary pixels
    fluor_filt=fluor_filt.*~boundaryRegion;
 
    %serach for offset that when added to centerline, maximizes the overlap
    %between the centerline and the fluor_filt image, starting with
    %xyOffset2
    [xyOffset2_min,~]=fminsearch(@(x) CLsearch(fluor_filt,CL2X,...
        CL2Y,show,x),xyOffset2);
    
    
    % do correlation to find xy offset between images
    % z project highmag and transfrome fluor images 
    hiResProj=normalizeRange(sum(worm_im,3));
    
    %do correlation to find offsets 
    corrIm=conv2(fluor_thresh_proj,rot90(hiResProj,2),'same');
    [CLoffsetY,CLoffsetX]=find(corrIm==max(corrIm(:)));
    CLoffsetX=CLoffsetX-round(size(fluor_thresh_proj,2)/2);
    CLoffsetY=CLoffsetY-round(size(fluor_thresh_proj,1)/2);  
    
    %add together all centerline offsets. 
    xyOffset3=xyOffset2_min-[CLoffsetX CLoffsetY];
    
    %apply offsets so that centerlines better overlap with high mag images
    CL2_hi_long(:,2,:)=CL2_hi_long(:,2,:)+xyOffset3(2);
    CL2_hi_long(:,1,:)=CL2_hi_long(:,1,:)+xyOffset3(1);
    
    if show
        close all
        imagesc(worm_im(:,:,midZ))
        hold on
        CL2X=CL2X+xyOffset3(1);
        CL2Y=CL2Y+xyOffset3(2);
        plot(CL2X,CL2Y,'xr');
    end
    
    %% crop cetnerline to retain region which is in the hi res image
    %find distance squared between centerline points in all frames and the
    %middle of the hi res image
    CLcenter=sum(bsxfun(@minus, CL2_hi_long,rect1(3:4)/2).^2,2);
    %find the index of the closest point. 
    [~,CLcenter]=min(CLcenter,[],1);
    CLcenter=mean(CLcenter(:));
    %add buffer to output length, will be cropped off later
    outputLength2=outputLength+outputLengthBuff;
    %add the output range to the closest point
    inImageRange=CLcenter+(-outputLength2:outputLength2);
    %interpolate to pull out CL points in that range.
    CL2_hi=interp1(CL2_hi_long,inImageRange,'*linear','extrap');
    CL2_size=size(CL2_hi);
    %% show middle imagee with overlayed centerline

    % pickout central centerline and middle slice
    % only needed for display purposes
    if show || isempty(side)
    %for middle image, take the 3 middle slices and zproject them.
    midIm=mean(worm_im(:,:,midZ+(-1:1)),3);
    midIm=double(midIm);
    
    %apply band pass filters
    midIm=bpass(midIm,2,[20,20]);
    
    %another badpass filter to really emphasize middle portion of
    %wormbrain. this helps pick out the central centerline
    midImS=imfilter(midIm,Sfilter);    
    
    %take the X and Y components of the centerlines
    CL3allX=(CL2_hi(:,1,:));
    CL3allY=(CL2_hi(:,2,:));
    
    %go through all centerlines and see which one has greates overlap with
    %middle slice midImS, that is the central centerline for visualization.
    minSearch=interp2(midImS,CL3allX,CL3allY);
    minSearch=squeeze(nansum(minSearch,1));
    minY=find(minSearch==max(minSearch));
    midZCL=round(mean(minY));
    CL_mid=CL2_hi(:,:,midZCL);
    end
    
    if show
        close all
        imagesc(midIm)
        hold on
        plot(CL_mid(:,1),CL_mid(:,2),'x')
        axis equal
        hold off
        drawnow
    end
    %% make coordinate system around the worm
    Tv=zeros(CL2_size(1),3,CL2_size(3));
    Bv=Tv; Nv=Tv;
    % for each centerline, make the Tangent, Normal, and Binormal vectors
    for iSlice=1:CL2_size(3);
        % T, N and B are first made in 2D,
        current_cl=CL2_hi(:,:,iSlice);
        T=normr(gradient(current_cl',5)');
        N=[T(:,2) -T(:,1)];
        B=T(:,1).*N(:,2)-T(:,2).*N(:,1);
        N=bsxfun(@times, N,sign(B));
        B=sign(B);
        % and then the Zcomponent is added
        T=[T zeros(size(current_cl(:,1)))];
        N=[N zeros(size(current_cl(:,1)))];
        B=[zeros(size(current_cl(:,1))) zeros(size(current_cl(:,1))) B];
        %compile results into vectors containing TBN for each centerline
        Tv(:,:,iSlice)=T;
        Nv(:,:,iSlice)=N;
        Bv(:,:,iSlice)=B;
        
    end
    %fix any sign flipping in the binormal and normal vector. Make the
    %binormal always be positive, tangent points from tail to head, normal
    %fits in according the right hand rule
    signVector=sign(Bv(:,3,:));
    Bv=bsxfun(@times,Bv,signVector);
    Nv=bsxfun(@times,Nv,signVector);
    
    %% select worm orientation 
    % the side input specifies the side of the nervcord. If the nerve chord
    % is on the left, flip the B and N vectors so that in interpolation, it
    % will come out on the right. This is so that different worms look
    % similar. If you don't care about flipping worms to look similar, this
    % part doesnt matter and any input, right or left, can be used and the
    % worms might flip (as is currently the case). Would be nice to do a
    % simple machine vision task of automatically finding the side in the
    % initializer. 
    
    if isempty(side)
        %if side is not input, use popup to pick side
        imagesc(midIm)
        choice = menu('Which Side is the nerve chord on?','Right','Left');
        if choice==2
            side='Left';
        else
            side='Right';
        end
    end
    %flip worm if side is left. 
    if ~isempty(strfind(side,'eft'))
        Bv=-Bv;
        Nv=-Nv;
    end

    
    plane_num=size(Tv,1);
    %make the first and last 'endround' tbn vectors the same so nothing
    %strange happens at the ends.
    
    %if number of points is larger than 10, use end round of 5. If smaller
    %than 10, use half the range, but something probably went wrong. 
    if plane_num>10
        endround=5;
    else
        endround=round(plane_num/2);
    end
    
    %replace TBN close to the ends, replacing them with the last one before
    %the endround.
    for iSlice=1:CL2_size(3);
        for i=1:endround
            Tv(plane_num-i+1,:,iSlice)=Tv(plane_num-endround+1,:,iSlice);
            Bv(plane_num-i+1,:,iSlice)=Bv(plane_num-endround+1,:,iSlice);
            Nv(plane_num-i+1,:,iSlice)=Nv(plane_num-endround+1,:,iSlice);
            Tv(i,:,iSlice)=Tv(endround,:,iSlice);
            Bv(i,:,iSlice)=Bv(endround,:,iSlice);
            Nv(i,:,iSlice)=Nv(endround,:,iSlice);
        end
    end
    

    
    
    %% show stack with centerline
    if show
        close all
        for iSlice=1:raw_size(3);
            imagesc(worm_im(:,:,iSlice));colormap hot
            hold on
            clSlice=iSlice;
            clSlice=round(clSlice);
            clSlice(clSlice<1)=1;
            clSlice(clSlice>CL2_size(3))=CL2_size(3);
            plot(CL2_hi(:,1,clSlice),CL2_hi(:,2,clSlice));
            quiver(CL2_hi(1:10:end,1,clSlice),CL2_hi(1:10:end,2,clSlice),...
                Nv(1:10:end,1,clSlice),Nv(1:10:end,2,clSlice))
            hold off
            axis auto equal off
            xlim([0 600]);ylim([0 600])
            pause(.1)
        end
    end
    
    %% build straightening coordinate system. 
    %define range to interpolate outward from centerline, adding a buffer
    %that will be removed later after additional template alignments. 
    

    
    outputRadius2=outputRadius+outputRadiusBuff;
    %create a 2*window +1 square around each point for interpolationg using
    %the B and N vectors
    [J,K]=meshgrid(-outputRadius2:outputRadius2,...
        -outputRadiusZ:outputRadiusZ);
    
    % build z coordinates using Nv and Bv vectors. For this, the binormal
    % vector Bv will almost always point in the Z direction, while the Z
    % component of the normal vector will be zero. 

    % the coordinate system is built by multiplying the J and K coordinates
    % with the normal and binormal vectors, taking only the Z direction.
    % This will tell me which slices in the unstraightened images to pull
    % from
    zslice=bsxfun(@times,J,permute(Nv(:,3,1),[3,2,1]))*zRatio+...
        bsxfun(@times,K,permute(Bv(:,3,1),[3,2,1]))*zRatio+midZ;
    
    % get the z values (which slice from the raw image) for every slice in 
    % the straightened space. They will all be the same so I just need one 
    % point from each slice. 
    zLevels=((zslice(:,1,1)));
    
    %correct for non monoticity, there is occasionally an error where the z
    %values are non monotonic because of small noise or are repeated.  
    if sign(nanmean(diff(zRange)))==-1
        % for downstroke of the triangle wave
        
        %cummin ensures that the range is never increasing, then take
        %unique so that the range is monotonic. 
        zRange2=unique(cummin(zRange)); 
        zRange2=zRange2-zRange2(midZ);
        %do interpolation trick to avoid any repeats
        adjusted_z=zLevels/10-zLevels(round(outputRadiusZ+1))/10;
        zInterp=interp1(zRange2,1:length(zRange2),adjusted_z);
        %flip z values so that stack reads from low to hi z.
        zInterp=flipud(zInterp);
        zLevels=flipud(zLevels);
    else
        %for the upstroke of triangle wave
        zRange2=unique(cummax(zRange));
        zRange2=zRange2-zRange2(midZ);
        %do interpolation trick to avoid any repeats

        adjusted_z=zLevels/10-zLevels(round(outputRadiusZ+1))/10;
        zInterp=interp1(zRange2,1:length(zRange2),adjusted_z);
    end
    
    % remake zslice with new zInterp so direction is correct with no
    % repeats
    zslice=repmat(zInterp,1,2*outputRadius2+1,size(Bv,1));
    
    %bound the z at the top by the image size and bottom by the lowest
    %unique centerline. 
    zLevels(zLevels<min(ia))=min(ia);
    zLevels(zLevels>raw_size(3))=raw_size(3);
    
    %interpolate the centerlines X and Y coordinates to get one for each
    %zLevel. 
    CL3xinterp=interp1(squeeze(CL2_hi(:,1,:))',zLevels,'linear')';
    CL3xinterp=permute(CL3xinterp,[2,3,1]);
    CL3yinterp=interp1(squeeze(CL2_hi(:,2,:))',zLevels,'linear')';
    CL3yinterp=permute(CL3yinterp,[2,3,1]);
    
    %interpolate the normal and binormal vector coordinate systems to get
    %one for each z level. 
    NvInterpx=interp1(squeeze(Nv(:,1,:))',zLevels,'linear');
    NvInterpy=interp1(squeeze(Nv(:,2,:))',zLevels,'linear');
    BvInterpx=interp1(squeeze(Bv(:,1,:))',zLevels,'linear');
    BvInterpy=interp1(squeeze(Bv(:,2,:))',zLevels,'linear');
    
    %build the coordinate systems by first muliplying the J and K space
    %along the N and V vectors, then adding that to the centerline.
    
    xslice=bsxfun(@times,J,permute(NvInterpx,[1,3,2]))+...
        bsxfun(@times,K,permute(BvInterpx,[1,3,2]));
    xslice=bsxfun(@plus, xslice,CL3xinterp);
    
    yslice=bsxfun(@times,J,permute(NvInterpy,[1,3,2]))+...
        bsxfun(@times,K,permute(BvInterpy,[1,3,2]));
    yslice=bsxfun(@plus, yslice,CL3yinterp);
    
    %permute transformations so sizes and directions work. 
    xslice=permute(xslice,[3,2,1]);
    yslice=permute(yslice,[3,2,1]);
    zslice=permute(zslice,[3,2,1]);
 
    % round all transformations
    xslice=round(xslice);yslice=round(yslice);zslice=round(zslice);
    
    %% Straighten worm image for correlation alignment with template
    % Straightening happens twice, this is the first time and is used to
    % get the transformations for image stabilization and alignment with
    % template. With those additional transformations, we will straighten
    % again later.
    
    %use only the pixels that are within the bounds of the original image
    inImageMap=xslice>0 & zslice>0 & yslice>0 & xslice<raw_size(2) &...
        yslice<raw_size(1) &  zslice<raw_size(3);
    inImageMapIdx=sub2ind_nocheck(raw_size,(yslice(inImageMap)),...
        (xslice(inImageMap)),(zslice(inImageMap)));
    
    % use indexing to apply the straightening transformation on the
    % unsmoothed images.
    V=zeros(size(xslice));
    V(inImageMap)=worm_im(inImageMapIdx);
    
    
    
    %% stack stabilization
    % get all transformations between slices, starting from bottom of stack
    % to the top. tformAll is a cell array with a transform between each
    % slice
    [~, tformAll]=stackStabilization(V,30,show,0);
    R = imref2d(size(V(:,:,1))) ;
    for iSlice=1:size(xslice,3)
        if any(any(inImageMap(:,:,iSlice)))
            %apply the transformations to both x and y transformations
            %slice by slice. Its faster to concatenate them, then call
            %imwarp on both of them.
            xyslice=cat(3,xslice(:,:,iSlice),yslice(:,:,iSlice));
            temp=imwarp(...
                xyslice, R,tformAll{iSlice},'nearest','OutputView',R);
            xslice(:,:,iSlice)=temp(:,:,1);
            yslice(:,:,iSlice)=temp(:,:,2);
        end
    end
    %round all transformations
    xslice=round(xslice);yslice=round(yslice);zslice=round(zslice);
    %use only the pixels that are within the bounds of the original image
    inImageMap=xslice>0 & zslice>0 & yslice>0 & xslice<raw_size(2) &...
        yslice<raw_size(1) &  zslice<raw_size(3);
    inImageMapIdx=sub2ind(raw_size,(yslice(inImageMap)),...
        (xslice(inImageMap)),(zslice(inImageMap)));
    
    Vsmooth=zeros(size(xslice));
    V=Vsmooth;
    % use indexing to apply the straightening transformation on the
    % smoothed and unsmoothed images. 
    Vsmooth(inImageMap)=worm_smooth(inImageMapIdx);
    V(inImageMap)=worm_im(inImageMapIdx);
    imsize=size(V);

    %% Correlation align with template image
    
    if ~isempty(Vtemplate)
        Vproj=squeeze(nansum(V,3)); %do z projection

        %do cross correlation and look for peak to find offset which
        %maximizes overlap between the template and the current Vproj
        xIm=conv2(Vproj,rot90(Vtemplate,2),'same');
        [xlag,ylag]=find(xIm==max(xIm(:)));
        lags=[xlag,ylag]-round(imsize(1:2)/2);
        
    else
        %if no template is given, assume no shifts are required.
        lags=[0 0];
    end
    %% apply offsets to the straightened image and all of the transformations
    %make ndgrid for interpolation
    [ndX,ndY,ndZ]=ndgrid(1:imsize(1),1:imsize(2),1:imsize(3));
    %apply offsets
    ndX=ndX+lags(1);
    ndY=ndY+lags(2);
    ndX=round(ndX);ndY=round(ndY);
    
    % get pixels that are still within the image after the shift that would
    % happen from the correlation alignment
    %something strange here, 
    inImage=(ndY>0 & ndX>0 &  ndX<=imsize(1) & ndY<=imsize(2));
    inImageIdx=sub2ind_nocheck(...
        imsize,ndX(inImage),ndY(inImage),ndZ(inImage));
    output_yrange=outputLengthBuff+1:imsize(1)-outputLengthBuff;
    output_xrange=outputRadiusBuff+1:imsize(2)-outputRadiusBuff;
    
    %use logical indexing (instead of interp) to shift the images
    temp=zeros(imsize);
    temp(inImage)=V(inImageIdx);
    %crop off the lenghth and radius buffers
    temp2=temp(output_yrange,output_xrange,:);
    V=temp2;
    % do the same for each of the transformations
    temp(inImage)=xslice(inImageIdx);
    temp2=temp(output_yrange,output_xrange,:);
    xslice=temp2;
    
    temp(inImage)=yslice(inImageIdx);
    temp2=temp(output_yrange,output_xrange,:);
    yslice=temp2;
    
    temp(inImage)=zslice(inImageIdx);
    temp2=temp(output_yrange,output_xrange,:);
    zslice=temp2;
    
    %do the same for the smoothed straightened image
    temp(inImage)=Vsmooth(inImageIdx);
    temp2=temp(output_yrange,output_xrange,:);
    Vsmooth=temp2;
    
    %% segmentation
    imsize=size(V);

    V(isnan(V))=0;
    Vsmooth(isnan(Vsmooth))=0; 
    % option, to use presmoothed version, much faster but may not be a s good
    % do segmentation on straightened worm
    [wormBW2,~]=WormSegmentHessian3dStraighten(V,options,Vsmooth);
    %% Filter out noise at top and bottom of stack
    % bpass filters tend to add noise at the top and the bottom of the
    % stack, we want to remove this. This will be done by projecting along
    % the x and y of the image mask to get a plot of number of pixels vs Z.
    % We expect to find a peak in the number of pixels close to the ends of
    % the stack. 
    
    %project along x and y
    BWplot=(squeeze(sum(sum(wormBW2,1),2)));
    %smooth
    BWplot=smooth(BWplot,20);
    %find peaks on  the projection plot
    [~,peak_loc]=findpeaks(BWplot);
    %find lowest and highest peaks, corresponding to the peaks in noise
    %near the edges
    endpts=peak_loc([1,length(peak_loc)]);
    %find the trough by finding peaks of -BWplot,
    [~,locs]=findpeaks(-BWplot);
    %find the ones closes the top and bottom of the worm, this will define
    %deltion zones between 1 and botpoint1, and botpoints2 and the end of
    %the stack
    botpoint1=locs((locs>endpts(1)));
    botpoint2=locs((locs<endpts(2)));
    % if no troughs are found, great, just use the top and bottom of the
    % stack (this will end up doing nothing in the end)
    if isempty(botpoint1);botpoint1=1;end;
    botpoint1=botpoint1(1);
    
    if isempty(botpoint2);botpoint2=imsize(3);end;
    botpoint2=botpoint2(end);
    % if the troughs founds are more than a quarter away from the edges,
    % then there probably isnt any noise artifact at the edges, so move 
    botpoint1(botpoint1>imsize(3)*1/4)=1;
    botpoint2(botpoint2<imsize(3)*3/4)=imsize(3);
    
    
    cc=bwconncomp(wormBW2,6);
    
    %kill of objects that touch the deletion zone
    %function that take linear pixel idx and returns Z
    zindexer=@(x,s) x./(s)+1; 
    %apply zindexer to each of the objects from bwconncomp
    objectZ=cellfun(@(x) zindexer(x,imsize(1)*imsize(2)),cc.PixelIdxList...
        ,'uniformoutput',0);
    
    %label objects as bad if theh touch the deletion zones
    badRegions_bot=cellfun(@(x) any(objectZ<=botpoint1),objectZ);
    badRegions_top=cellfun(@(x) any(objectZ>=botpoint2),objectZ);
    badRegions=(badRegions_bot | badRegions_top)';
    
    % remove bad objects by setting all their pixels to false
    wormBW2(cell2mat(cc.PixelIdxList(badRegions)'))=false;
    cc.PixelIdxList=cc.PixelIdxList(~badRegions);
    cc.NumObjects=nnz(~badRegions);
    
    
    %% compile results into pointStats structure for saving
    %hard cap at 200 neurons, occasionally  you have big fails that produce
    %many hundreds of points. this is bad, just blank everything if this
    %happens
    cc=bwconncomp(wormBW2,6);

    if cc.NumObjects<200
        %get results from binary mask using region props
        stats=regionprops(cc,V,'Centroid','MeanIntensity',...
            'Area');
        %get roi mean intensities
        intensities=[stats.MeanIntensity]';
        %get roi centroids
        P=cell2mat({stats.Centroid}');
        
        %flip x and y
        P(:,[1 2])=P(:,[2 1]);
        
        %get roi volumes
        Areas=[stats.Area]';
        %centroids from unstraightened coordinate system using
        %interpolation. 
        Poriginal=[interp3(xslice,P(:,2),P(:,1),P(:,3)) ...
            interp3(yslice,P(:,2),P(:,1),P(:,3)) ...
            interp3(zslice,P(:,2),P(:,1),P(:,3))];
        
        % load up all of these into pointstats
        pointStats.straightPoints=P(:,1:3);
        pointStats.rawPoints=Poriginal;
        pointStats.pointIdx=(1:cc.NumObjects)';
        pointStats.Rintensities=intensities;
        pointStats.Volume=Areas;
        
    else
        % if more than 200 points, make all of these empty
        pointStats.straightPoints=[];
        pointStats.pointIdx=[];
        pointStats.Rintensities=[];
        pointStats.Volume=[];
        pointStats.rawPoints=[];
    end
    % save stack index
    pointStats.stackIdx=iStack;
    % save the binary mask
    pointStats.baseImg=logical(wormBW2);
    % save coordinate lookup tables for going back and to unstraightened
    % coordinate system
    pointStats.transformx=uint16(xslice);
    pointStats.transformy=uint16(yslice);
    pointStats.transformz=uint16(zslice);

    
    %% save results
    %set up paths for saving
    image_destination=...
        [destination_path filesep 'image' num2str(iStack,'%3.5d') '.tif'];
    mat_destination=...
        [destination_path filesep 'pointStats' num2str(iStack,'%3.5d')];
    %if show is set to 1, save aditional centerlines for and volume for
    %easy visualization of straightening later. Never happens on cluster
    if show>1
        save_destination=...
            [destination_path filesep 'saveFile' num2str(iStack,'%3.5d')];
        save(save_destination,'CL2_hi_long', 'CL2_hi','Tv','Bv','Nv','worm_im',...
            'V','Vsmooth','pointStats');
    end
    %save tif
    tiffwrite(image_destination,single(V),'tif');
    % save pointStats into mat file
    save(mat_destination,'pointStats');

    fclose(Fid);
catch me
    %if an error occurs, display the error and save the entire workspace
    %for future trouble shooting
    display(me.identifier)
    err_destination=...
        [destination_path filesep 'ERROR' num2str(iStack,'%3.5d')];
    save(err_destination,'me');
end
