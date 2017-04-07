
function fiducialCropper3( dataFolder)
%%%%%%%%%%%
% Fiducialcropper3 takes dataFolder as input after it has gone through all
% of the wormAnalysisPipeline. It takes the pointStats file with all of the
% coordinates of the tracked neurons in the straightened coordinates and
% extracts the Red and Green signals from the original .dat file. It also
% analyzes the behavior images and extracts the behavior. Outputs are saved
% into a heatData.mat file in the dataFolder

%% create some Gaussian Kernals


filterKernal2=makeGaussKernel(500);


[RvalAll,GvalAll]=extractSignal(dataFolder);

%% makes photobleaching corrected, spatially corrected heatmaps, saves result
heatMapGeneration(dataFolder,RvalAll,GvalAll)
load([dataFolder filesep 'heatData'])

%% make timing tracks

[bfAll,~,hiResData]=tripleFlashAlign(dataFolder);

hasPoints=1:length(pointStats);
hasPointsTime=hiResData.frameTime(diff(hiResData.stackIdx)==1);
hasPointsTime=hasPointsTime(hasPoints);


[centerline,eigenProj, CLV,wormCentered]=loadCLBehavior(dataFolder);
n_CL=size(centerline,3);
clTime=bfAll.frameTime(1:n_CL);



%% load centerline data, pad with zeros if needed to making size the same as
%behavior video size

CLbehavior=sign(smooth(-CLV',50));
hiResCLbehavior=interp1(clTime,CLbehavior,hiResData.frameTime,'nearest',0);

%% add centerline center of mass to stage position
[xPos,yPos]=calculate_cm_position(centerline,hiResData,clTime);
%get worm cm positions for each volume measured


%% Calculate worm velocities

filterFactor2=imfilter(ones(size(xPos)),filterKernal2);
xV=imfilter((xPos),filterKernal2)./filterFactor2;
yV=imfilter((yPos),filterKernal2)./filterFactor2;
hiResVel=[gradient(xV) gradient(yV)];
hiResVel=sqrt(sum(hiResVel.^2,2));

hiResVel=hiResVel.*hiResCLbehavior;
hiResVel=smooth(hiResVel,500);
%%
xPosTrack=xPos((find(diff(hiResData.stackIdx)>0)));
yPosTrack=yPos((find(diff(hiResData.stackIdx)>0)));
xPosTrack=xPosTrack(hasPoints);
yPosTrack=yPosTrack(hasPoints);

vTrack=hiResVel((find(diff(hiResData.stackIdx)>0)));
vTrack=vTrack(hasPoints);
%% load eigen behavior and do cylindrical rotation for more independent Z

eigenBehavior=eigenProj;
temp=eigenBehavior;
behaviorProjection=temp(1:3,:)';
[~,maxIdx]=max(sum(behaviorProjection.^2,2));
maxPoint=behaviorProjection(maxIdx,:);
[phi,theta,R] = cart2sph(maxPoint(1),maxPoint(2),maxPoint(3));
[x,E]=fminsearch(@(x) circlePlane(behaviorProjection,x(1),x(2)),[phi theta]);
phi=x(2);
theta=x(1);
[~,projections]=circlePlane(behaviorProjection,theta,phi);
n=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
behaviorZ=behaviorProjection*n';
Rmat=normc(projections'*behaviorZ);
Rmat=[Rmat(1) Rmat(2); Rmat(2) -Rmat(1)];
projections=(Rmat*projections')';


behaviorZTrack=interp1(clTime,behaviorZ,hasPointsTime);
projectionsTrack=interp1(clTime,projections,hasPointsTime);
hiResBehaviorZ=interp1(clTime,behaviorZ,hiResData.frameTime);


%% calculate angles in PC1 and PC2 space
behaviorTheta=atan2(projections(:,1),projections(:,2));
behaviorTheta=unwrap(behaviorTheta);
behaviorTheta=imfilter(behaviorTheta,filterKernal);



%% make some sort of ethogram;
%still in progress for more types of worm
% -1 for reverse
% 0 for pause
% 1 for forward
% 2 for turn

ethogram=makeEthogram(hiResVel,hiResBehaviorZ);
ethoTrack=interp1(hiResData.frameTime,ethogram,hasPointsTime,'nearest');

%combine behaviors into structure
behavior.ethogram=ethoTrack;
behavior.x_pos=xPosTrack;
behavior.y_pos=yPosTrack;
behavior.v=vTrack;
behavior.pc1_2=projectionTrack;
behavior.pc_3=behaviorZTrack;

%% save files
save([dataFolder filesep 'positionData'],...
    'xPos','yPos','hiResVel','behaviorZTrack','projectionsTrack')

save([dataFolder filesep 'heatData'],...
    'behavior',...
    'hasPointsTime',...
    'trackWeightAll',...
    'bfTime',...
    '-append');



function boxPix=boxCoord(L)
%% initialize box. This will be added to gather pixel values from around the centroid
boxR=[L L L];
[boxX,boxY,boxZ]=meshgrid(-boxR(1):boxR(1),-boxR(2):boxR(2),-boxR(3):boxR(3));
boxMask=(boxX.^2+boxY.^2+boxZ.^2)<25;
boxX=boxX(boxMask);
boxY=boxY(boxMask);
boxZ=boxZ(boxMask);
boxPix=[boxX,boxY,boxZ];


function all_corr=illuminationCorrection()
%% create illumination profile correction, if files are present

try
    %load illumination profiles, these are images normalized images of a
    %fluorescent calibration slide. Each image is full field, but due to
    %the filters only the appropriate half of the dual view image is shown.
    
    profileG=load('illumination_profile_G.mat');
    profileG=profileG.illumination_profile;
    g_corr=1./profileG;
    %remove very bright and very dim pixels in the calibration image
    g_corr(g_corr>5| g_corr<0)=0;
    
    profileR=load('illumination_profile_R.mat');
    profileR=profileR.illumination_profile;
    r_corr=1./profileR;
    r_corr(r_corr>5 | r_corr<0)=0;
    
    %combine the two halves, the two regions corresponding to the image should
    %have no overlap so a straight pix by pix sum will work
    all_corr=g_corr+r_corr;
catch me
    
    display(' No illumination profile found, no correction applied')
    all_corr=1;
end


function [RvalAll,GvalAll]=extractSignal(dataFolder)
%%
%get files from dataFolder, including the imageFolder and the pointStats
%file.

boxPix=boxCoord(5);
all_corr=illuminationCorrection();

pointStatsFile=[dataFolder filesep 'pointStatsNew.mat'];
psFolder=dir([dataFolder filesep 'CLstraight*']);
imageFolder=[dataFolder filesep psFolder(end).name];
%% load pointStats File

pointStats=load(pointStatsFile);
pointStats=pointStats.pointStatsNew;

% read in sCMOS dat file
datFileDir=dir([dataFolder filesep 'sCMOS_Frames_U16_*']);
datFile=[dataFolder filesep datFileDir.name];
%get image dimensions
[rows,cols]=getdatdimensions(datFile);
nPix=rows*cols;

%read in timing data
[~,~,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);


%also need phase delay used between z position and image number
timeOffset=load([dataFolder filesep 'startWorkspace.mat'],'zOffset');
timeOffset=timeOffset.zOffset;

% load alignments
alignments=load([dataFolder filesep 'alignments.mat']);
alignments=alignments.alignments;
S2AHiRes=alignments.S2AHiRes;
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
background=alignments.background;



%%
fiducialPoints=cell(1,length(pointStats));
for i=1:length(pointStats)
    try
        %load and process the image
        
        iFile=pointStats(i).stackIdx;
        Fid=fopen(datFile);
        %get startframe for single voluem scan
        lowLimit=find(hiResData.stackIdx==iFile,1,'first')+timeOffset;
        nSlices=nnz(hiResData.stackIdx==iFile);
        
        %move pointer and read the pixel values
        status=fseek(Fid,2*(lowLimit)*nPix,-1);
        pixelValues=fread(Fid,nPix*nSlices,'uint16',0,'l');
        
        %create image volume
        hiResImage=(reshape(pixelValues,rows,cols,nSlices));
        
        %subtract background and apply intensity correction
        hiResImage=bsxfun(@minus,hiResImage,background);
        hiResImage=bsxfun(@times,hiResImage,all_corr);
        hiResImage(hiResImage<0)=0;
        
        
        %% load transformation from striaghtened to original
        pointFile=([imageFolder ...
            filesep 'pointStats' num2str(iFile,'%2.5d') '.mat']);
        pointStatsTemp=load(pointFile);
        pointStatsTemp=pointStatsTemp.pointStats;
        
        %tgese are the lookuptables between pixel in straightened
        %coordinate and original image.
        
        X=double(pointStatsTemp.transformx);
        Y=double(pointStatsTemp.transformy);
        Z=double(pointStatsTemp.transformz);
        
        if any(size(X)==1);
            [X, Y, Z]=meshgrid(X,Y,Z);
        end
        
        nSize=size(X);
        
        %pull tracking results
        trackIdx=pointStats(i).trackIdx;
        trackWeights=pointStats(i).trackWeights;
        straightPoints=pointStats(i).straightPoints;
        
        % get the tracked neurons present in this volume
        present=~isnan(trackIdx) & ~isnan(straightPoints(1:length(trackIdx),1));
        
        %compress and organize tracking results
        trackIdx=trackIdx(present);
        trackWeights=trackWeights(present(1:length(trackWeights)));
        trackWeightstemp=zeros(size(trackIdx));
        trackWeightstemp(1:length(trackWeights))=trackWeights;
        trackWeights=trackWeightstemp;
        straightPoints=straightPoints(present,:);
        
        %% get pixels of interest by adding the boxPix to each point.
        boxPixAll=bsxfun(@plus,straightPoints,permute(boxPix,[3,2,1]));
        boxPixAll=round(boxPixAll);
        % bound at the end of the box, though this is rare
        boxPixAll(boxPixAll<1)=1;
        for iDim=1:3
            temp=boxPixAll(:,iDim,:);
            temp(temp>nSize(iDim))=nSize(iDim);
            boxPixAll(:,iDim,:)=temp;
        end
        
        %turn sub into linear index,
        boxPixAllLin=sub2ind_nocheck(nSize,boxPixAll(:,1,:),boxPixAll(:,2,:),boxPixAll(:,3,:));
        boxPixAllLin=permute(boxPixAllLin,[3,1,2]);
        %%
        %transform points ball points into unstraightened coord (red image)
        newX=X(boxPixAllLin);
        newY=Y(boxPixAllLin);
        newZ=Z(boxPixAllLin);
        newX1=newX+rect1(1)-1;
        newY1=newY+rect1(2)-1;
        
        %do transformation for green image
        [newX2,newY2]=transformPointsInverse(S2AHiRes.t_concord,newX,newY);
        newY2=newY2+rect2(2)-1;
        newX2=newX2+rect2(1)-1;
        newX2=round(newX2);
        newY2=round(newY2);
        %make raw points for saving, need to check
        rawPoints=coordinateTransform3d(straightPoints,X,Y,Z);
        hiResRange=find(hiResData.stackIdx==pointStats(i).stackIdx);
        hiResIdx=hiResRange(rawPoints(:,3))+timeOffset;
        hiResVoltage=interp1(hiResData.Z,hiResIdx-timeOffset);
        fiducialPointsi=cell(200,4);
        %% loop over points and get average pixel intensities
        %intialize intensity vector for a given volume
        
        R_i=interpSignal(hiResImage,newX1,newY1,newZ);
        G_i=interpSignal(hiResImage,newX2,newY2,newZ);
        
        
        %on first pass, initialize outputs for green and red and weights as
        %nans
        if ~exist('GvalAll','var')
            GvalAll=nan(nnz(present),max([pointStats.stackIdx]));
            RvalAll=GvalAll;
            trackWeightAll=RvalAll;
        end
        
        if length(trackIdx)==length(trackWeights)
            trackWeightAll(trackIdx,iFile)=trackWeights;
        end
        
        %save intensity averages
        RvalAll(trackIdx,iFile)=R_i;
        GvalAll(trackIdx,iFile)=G_i;
        
        %save cell aray with all coordinates, voltages, and frame numbers
        fiducialPointsi(trackIdx,:)=...
            num2cell([rawPoints(:,1:2) hiResVoltage hiResIdx]);
        
        fclose(Fid);
        fiducialPoints{iFile}=fiducialPointsi;
        display(['Finished frame : ' num2str(i)]);
    catch ME
        ME
        display(['Error frame:' num2str(i)]);
    end
    
end

%% save cropped regions and a new set of unstraighted centroids
newFiducialFolder=[dataFolder filesep 'BotfFiducialPoints'];
mkdir(newFiducialFolder);
%save time offset and unstraightened fiducial cell structure, not as
%important any more, but good for visualization
clickPoints=0;
save([newFiducialFolder filesep 'timeOffset'],'timeOffset');
save([newFiducialFolder filesep 'botFiducials'],...
    'fiducialPoints',...
    'clickPoints');


% interpSignal takes coordinates given by XYZ and interpolates the signal
% from hiResImage. The signal is averaged over the columns, with a few
% checks for ensuring the pixels are in the image and do not have wild
% variations.
function signal=interpSignal(hiResImage,X,Y,Z)

signal=nan(size(X,2),1);


for iPoint=1:size(X,2)
    %%
    %get all pixels of interest for a point iPoint, red version
    xball1=X(:,iPoint);
    yball1=Y(:,iPoint);
    zball=Z(:,iPoint);
    
    %remove bad pixels
    xball1(xball1==0)=nan;
    yball1(yball1==0)=nan;
    zball(zball==0)=nan;
    reject1=isnan(xball1) |isnan(yball1)|isnan(zball);
    keep=~reject1;
    
    % get index for each pixel of interest for R and G
    linBall1=sub2ind_nocheck(...
        size(hiResImage),...
        yball1(keep),...
        xball1(keep),...
        zball(keep));
    
    % get pixel values
    Fset=hiResImage(linBall1);
    %some pixels are prone to producing large values, try to remove
    %these
    bad_pix=abs(zscore(Fset))>5;
    Fset(bad_pix)=[];
    Fnorm=normalizeRange(Fset);
    Fthresh=Fnorm<graythresh(Fnorm);
    Fset(Fthresh)=[];
    %take trim mean to try to clean signal a bit
    F_point=trimmean(Fset,20);
    
    signal(iPoint)=F_point;
    
end


function [xPos,yPos]=calculate_cm_position(centerline,hiResData,bf_frameTime)


%some conversion factors
pos2mm=1/10000;
CL2mm=1/557; % 1mm per 557 pixels
%angle between stage positions and behavior camera
stageCamAngle=90;
stageCamAngle=stageCamAngle*pi/180;


%get mean of CL coordinates for center of mass
centerLinePosition=squeeze(mean(centerline,1));
centerLinePosition(centerLinePosition==0)=nan;
centerLinePosition=inpaint_nans(centerLinePosition);
centerLinePosition=bsxfun(@minus,centerLinePosition,mean(centerLinePosition,2));

%interpolate into imaging rate of hiRes
hiResCLposition=...
    [interp1(bf_frameTime,centerLinePosition(:,1),...
    hiResData.frameTime,'nearest',0),...
    interp1(bf_frameTime,centerLinePosition(:,2),...
    hiResData.frameTime,'nearest',0)];

hiResCLposition(isnan(hiResCLposition(:,1)),1)=nanmean(hiResCLposition(:,1));
hiResCLposition(isnan(hiResCLposition(:,2)),2)=nanmean(hiResCLposition(:,2));

hiResCLposition=hiResCLposition*CL2mm;


%rotation matrix that takes motion of stage direction to motion in low mag
%behavior image
Rmatrix=[-cos(stageCamAngle) sin(stageCamAngle);...
    -sin(stageCamAngle) -cos(stageCamAngle)];
hiResCLposition=(Rmatrix'*hiResCLposition')';

%add CL center of mass to stage coordinate
xPosStage=hiResData.xPos*pos2mm;
xPosStage([0; diff(xPosStage)]==0)=nan;
xPosStage=inpaint_nans(xPosStage);
yPosStage=hiResData.yPos*pos2mm;
yPosStage([0; diff(yPosStage)]==0)=nan;
yPosStage=inpaint_nans(yPosStage);

% switched some signs 1 and 2 for jeff cls, may need switching for some
% camera rotations
xPos=xPosStage-1*hiResCLposition(:,1);
yPos=yPosStage+1*hiResCLposition(:,2);

%filter positions
filterKernal=makeGaussKernel(50);
filterFactor=imfilter(ones(size(xPos)),filterKernal);

xPos=imfilter(xPos,filterKernal)./filterFactor;
yPos=imfilter(yPos,filterKernal)./filterFactor;



%simple function to make normalized gaussian kernel with sdev
function filterKernel=makeGaussKernel(sdev)
filterKernel=gausswin(sdev);
filterKernel=filterKernel-min(filterKernel(:));
filterKernel=filterKernel/sum(filterKernel);




% -1 for reverse
% 0 for pause
% 1 for forward
% 2 for turn

%colors for plotting
function c_out=ethocolormap()
fcolor=[0 1 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
turncolor=[0 0 1];%[117 112 179]/256;
pausecolor=[255 217 50]/256;
c_out=[bcolor;pausecolor;fcolor;turncolor];


function ethogram=makeEthogram(hiResVel,hiResBehaviorZ)
%attempt to make an ethogram of behavior using the velocity of the animal
%(which will have a sign to take into account forward and reversal) and the
%projection onto the "3rd" eigenworm, where peaks correspond to high
%curvature and likely turns. 

%use 3rd eigenmode for finding turns
ethoZ=hiResBehaviorZ;
ethoZ=ethoZ-nanmean(ethoZ);
ethoZ(isnan(ethoZ))=0;

%cluster velocities to find forward and back
Vcluster=smooth(hiResVel,100);
idx=sign(Vcluster);
idx(abs(Vcluster)<.00005)=0;
ethogram=idx;
ethogram((((abs(ethoZ)>2*std(ethoZ))& (ethogram>=0)))|abs(ethoZ)>10)=2;

% Kill off short behaviors unless they are reversals
for iBehavior=[0 1 2];
    cpause=bwconncomp(ethogram==iBehavior);
    shortPause=cellfun(@(x) length(x), cpause.PixelIdxList);
    shortPause= shortPause'<500;
    shortPause=cpause.PixelIdxList(shortPause);
    shortPause=cell2mat(shortPause');
    ethogram(shortPause)=nan;
    ethogram=colNanFill(ethogram,'nearest');
end

ethogram(hiResVel==0)=nan;


