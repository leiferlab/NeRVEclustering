%%%%%%%%%%%
% Fiducial cropper 2 takes straightened coordinates form a pointStats file
% and crops out regions corresponding to the unstraighted form. It then
% averages over those ROIs to get the intensity at that time. 
function fiducialCropper3( dataFolder)
if nargin==0
display('Select PointStatsFile');
pointStatsFile=uipickfiles();
pointStatsFile=pointStatsFile{1};
dataFolder=fileparts(pointStatsFile);
display('Select Image Folder');
imageFolder=uipickfiles('filterspec',dataFolder);
imageFolder=imageFolder{1};

else
    pointStatsFile=[dataFolder filesep 'pointStatsNew.mat'];
    psFolder=dir([dataFolder filesep 'CLstraight*']);
    imageFolder=[dataFolder filesep psFolder.name];
end
load(pointStatsFile)


filterKernal=gausswin(50);
filterKernal=filterKernal-min(filterKernal(:));
filterKernal=filterKernal/sum(filterKernal);

filterKernal2=gausswin(500);
filterKernal2=filterKernal2-min(filterKernal2(:));
filterKernal2=filterKernal2/sum(filterKernal2);
filterFactor2=imfilter(ones(size(xPos)),filterKernal2);

%%
pointStats2=pointStatsNew;

presentIdx=1:length(pointStats2);
rows=1200;
cols=600;
nPix=rows*cols;
%%
try
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);
catch
    if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));
    end
end

zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
timeOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));


%%
%%
try
    alignments=load([dataFolder filesep 'alignments.mat']);
    alignments=alignments.alignments;
    S2AHiRes=alignments.S2AHiRes;
    rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
catch
display('Select Hi Res Alignment')

S2AHiResFile=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiResFile{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
end

%%
boxR=[5 5 5];
[boxX,boxY,boxZ]=meshgrid(-boxR(1):boxR(1),-boxR(2):boxR(2),-boxR(3):boxR(3));
boxMask=(boxX.^2+boxY.^2+boxZ.^2)<25;
boxX=boxX(boxMask);
boxY=boxY(boxMask);
boxZ=boxZ(boxMask);
boxPix=[boxX,boxY,boxZ];



%%
fiducialPoints=cell(1,length(pointStats2));
for i=1:length(pointStats2)
   %      pointStats2(presentIdx(i)).transitionMatrix=TrackMatrix{i};
   
    try
        %load the image
            iFile=pointStats2(i).stackIdx;
        pointFile=([imageFolder  filesep 'pointStats' num2str(iFile,'%2.5d') '.mat']);

            Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);
        Rout=zeros(12,1);
        Gout=Rout;
        backgroundOut=Gout;
       
    lowLimit=find(hiResData.stackIdx==iFile,1,'first')+timeOffset;
     nSlices=nnz(hiResData.stackIdx==iFile);
   
    status=fseek(Fid,2*(lowLimit)*nPix,-1);
    pixelValues=fread(Fid,nPix*nSlices,'uint16',0,'l');
    
    

    hiResImage=(reshape(pixelValues,rows,cols,nSlices));
    
    
        %% load transformation from striaghtened to original
    
     pointStatsTemp=load(pointFile);
pointStatsTemp=pointStatsTemp.pointStats;
X=double(pointStatsTemp.transformx);
Y=double(pointStatsTemp.transformy);
Z=double(pointStatsTemp.transformz);
    if any(size(X)==1);
        [X, Y, Z]=meshgrid(X,Y,Z);
    end
    nSize=size(X);
    
trackIdx=pointStats2(i).trackIdx;
    straightPoints=pointStats2(presentIdx(i)).straightPoints;
  %  rawPoints=rawPoints(~isnan(pointStats2(presentIdx(i)).trackIdx),:);
    present=~isnan(trackIdx) & ~isnan(straightPoints(1:length(trackIdx),1));
    trackIdx=trackIdx(present);
    if i==1
        GvalAll=nan(nnz(present),max([pointStats2.stackIdx]));
        RvalAll=GvalAll;
    end

    straightPoints=straightPoints(present,:);
    boxPixAll=bsxfun(@plus,straightPoints,permute(boxPix,[3,2,1]));
    boxPixAll=round(boxPixAll);
    boxPixAll(boxPixAll<1)=1;
    for iDim=1:3
        temp=boxPixAll(:,iDim,:);
        temp(temp>nSize(iDim))=nSize(iDim);
        boxPixAll(:,iDim,:)=temp;
    end
    
    boxPixAllLin=sub2ind_nocheck(nSize,boxPixAll(:,1,:),boxPixAll(:,2,:),boxPixAll(:,3,:));
    boxPixAllLin=permute(boxPixAllLin,[3,1,2]);
    %%
    %transform points ball points into unstraightened coord
    newX=X(boxPixAllLin);
    newY=Y(boxPixAllLin);
    newZ=Z(boxPixAllLin);
    newX1=newX+rect1(1)-1;
    newY1=newY+rect1(2)-1;    
   [newX2,newY2]=transformPointsInverse(S2AHiRes.t_concord,newX,newY);
    newY2=newY2+rect2(2)-1;
    newX2=newX2+rect2(1)-1;

    rawPoints=coordinateTransform3d(straightPoints,X,Y,Z);
    rawPoints(rawPoints<1)=1;
    hiResRange=find(hiResData.stackIdx==pointStats2(i).stackIdx);
    hiResIdx=hiResRange(rawPoints(:,3))+timeOffset;
    hiResV=interp1(hiResData.Z,hiResIdx-timeOffset);
    fiducialPointsi=cell(200,4);
    
    newHiIdx=interp1(hiResRange,newZ,'*nearest')+timeOffset;
    newV=interp1(hiResData.Z,newHiIdx-timeOffset);
    %%
    Rtemp=nan(length(trackIdx),1);
    Gtemp=nan(length(trackIdx),1);
    
    
    for iPoint=1:size(newX,2)
       %%
        xball1=newX1(:,iPoint);
        yball1=newY1(:,iPoint);
        zball=newZ(:,iPoint);
        xball1(xball1==0)=nan;
        yball1(yball1==0)=nan;
        zball(zball==0)=nan;
        reject1=isnan(xball1) |isnan(yball1)|isnan(zball);
        keep1=~reject1;
        
        xball2=round(newX2(:,iPoint));
        yball2=round(newY2(:,iPoint));
        xball2(xball2==0)=nan;
        yball2(yball2==0)=nan;
        reject2=isnan(xball2) |isnan(yball2)|isnan(zball);
        keep2=~reject2;
        
        
        linBall1=sub2ind_nocheck(size(hiResImage),yball1(keep1),xball1(keep1),zball(keep1));
        linBall2=sub2ind_nocheck(size(hiResImage),yball2(keep2),xball2(keep2),zball(keep2));
        Rset=hiResImage(linBall1);
        Gset=hiResImage(linBall2);
        Rset(Rset<100)=[];
        Gset(Gset<100)=[];
        Rnorm=normalizeRange(Rset);
        Rset(Rnorm<graythresh(Rnorm))=[];
        Gnorm=normalizeRange(Gset);
        Gset(Gnorm<graythresh(Gnorm))=[];
        Rval=trimmean(Rset,20);
        Gval=trimmean(Gset,20);
        
        Rtemp(iPoint)=Rval;
        Gtemp(iPoint)=Gval;
        
    end
    
        RvalAll(trackIdx,iFile)=Rtemp;
        GvalAll(trackIdx,iFile)=Gtemp;
         %   rawPointsCell=num2cell(rawPoints);
     fiducialPointsi(trackIdx,:)=num2cell([rawPoints(:,1:2) hiResV hiResIdx]);

fclose(Fid);
     fiducialPoints{iFile}=fiducialPointsi;
    display(['Finished frame : ' num2str(presentIdx(i))]);
    catch ME
        ME
        display(['Error frame:' num2str(presentIdx(i))]);
    end
    
end


%% save cropped regions and a new set of unstraighted centroids
newFiducialFolder=[dataFolder filesep 'BotfFiducialPoints'];
mkdir(newFiducialFolder);
clickPoints=0;
save([newFiducialFolder filesep 'timeOffset'],'timeOffset');

save([newFiducialFolder filesep 'botFiducials'],'fiducialPoints','clickPoints');
save([dataFolder filesep  'WormFiducialIntensityPoints_thresh'],'RvalAll','GvalAll')

%% makes photobleaching corrected, spatially corrected heatmaps, saves result
heatMapGeneration_spaceCorr(dataFolder,fiducialPoints,alignments,[])
load([dataFolder filesep 'heatData0'])
[bfTime,hiResFrameTime,hasPoints,bfRange,hiResRange,hasPointsTime,lookup]=dataTimeAlignment(dataFolder,fiducialPoints);
[centerline, offset,eigenProj, CLV,wormCentered]=loadCLBehavior(dataFolder);

frameTimeTrack=hiResData.frameTime((find(diff(hiResData.stackIdx)>0)));
frameTimeTrack=frameTimeTrack(hasPoints);
firstTime=nanmin(frameTimeTrack);
frameTimeTrack=frameTimeTrack-firstTime;


%% load centerline data
if offset>0
v= [zeros(1,offset)  CLV]';
else
    v=CLV';
end
zeroPad=zeros(size(bfAll.frameTime,1)-length(v),1);
v=[v;zeroPad];
v(1)=0;
CLV=v;
CLbehavior=sign(smooth(-CLV,50));
hiResCLbehavior=interp1(bfAll.frameTime,CLbehavior,hiResData.frameTime,'nearest','extrap');
hiResCLV=interp1(bfAll.frameTime,CLV,hiResData.frameTime,'nearest','extrap');

%% add centerline center of mass to stage position 
centerLinePosition=squeeze(mean(centerline,1));
centerLinePosition(centerLinePosition==0)=nan;
centerLinePosition=inpaint_nans(centerLinePosition);
centerLinePosition=bsxfun(@minus,centerLinePosition,mean(centerLinePosition,2));
v=zeros(2,length(CLV));
v(:,offset+1:offset+length(centerLinePosition))=centerLinePosition;
zeroPad=zeros(size(bfAll.frameTime,1)-length(v),1);
v=[v;zeroPad];
centerLinePosition=v';
hiResCLposition=...
    [interp1(1:length(centerLinePosition),centerLinePosition(:,1),...
        lookup.BF2hiRes,'nearest','extrap'),...
    interp1(1:length(centerLinePosition),centerLinePosition(:,2),...
        lookup.BF2hiRes,'nearest','extrap')];

hiResCLposition(isnan(hiResCLposition(:,1)),1)=nanmean(hiResCLposition(:,1));
hiResCLposition(isnan(hiResCLposition(:,2)),2)=nanmean(hiResCLposition(:,2));

hiResCLposition=hiResCLposition*1/557; % 1mm per 557 pixels
stageCamAngle=90;
stageCamAngle=stageCamAngle*pi/180;
%rotation matrix that takes motion of stage direction to motion in low mag
%behavior image
Rmatrix=[-cos(stageCamAngle) sin(stageCamAngle);...
    -sin(stageCamAngle) -cos(stageCamAngle)];
hiResCLposition=(Rmatrix'*hiResCLposition')';

xPosStage=hiResData.xPos/10000;
xPosStage([0; diff(xPosStage)]==0)=nan;
xPosStage=inpaint_nans(xPosStage);
yPosStage=hiResData.yPos/10000;
yPosStage([0; diff(yPosStage)]==0)=nan;
yPosStage=inpaint_nans(yPosStage);
xPos=xPosStage-1*hiResCLposition(:,1); % switched some signs 1 and 2 for jeff cls
yPos=yPosStage+1*hiResCLposition(:,2);
yPosStageTrack=yPosStage(hiResRange);
xPosStageTrack=xPosStage(hiResRange);
yPosStageTrack=yPosStageTrack(hasPoints);
xPosStageTrack=xPosStageTrack(hasPoints);
filterFactor=imfilter(ones(size(xPos)),filterKernal);
xV=imfilter(xPos,diff(filterKernal));
 xPos=imfilter(xPos,filterKernal)./filterFactor;
 yPos=imfilter(yPos,filterKernal)./filterFactor;

hiResxPos=xPos(hiResRange);
hiResyPos=yPos(hiResRange);
%center coordinates around zero
hiResxPos=hiResxPos-(max(hiResxPos)+min(hiResxPos))/2;
hiResyPos=hiResyPos-(max(hiResyPos)+min(hiResyPos))/2;

%% Calculate worm velocities

xPosTrack=xPos((find(diff(hiResData.stackIdx)>0)));
xPosTrack=xPosTrack(hasPoints);
yPosTrack=yPos((find(diff(hiResData.stackIdx)>0)));
yPosTrack=yPosTrack(hasPoints);
 xV=imfilter((xPos),filterKernal2)./filterFactor2;
 yV=imfilter((yPos),filterKernal2)./filterFactor2;
v=[gradient(xV) gradient(yV)];
v=sqrt(sum(v.^2,2));
D=[cumsum(v)];

hiResV=v;%gradient(D,50)./median(gradient(hiResData.frameTime));
hiResV=hiResV.*hiResCLbehavior;
hiResV=smooth(hiResV,500);
%% load eigen behavior and do cylindrical rotation for more independent Z

eigenBehavior=zeros(size(eigenProj,1),length(CLV));
eigenBehavior(:,offset+1:offset+length(eigenProj))=eigenProj;
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
hiResBehaviorZ=interp1(bfAll.frameTime,behaviorZ,hiResData.frameTime);


%% calculate angles in PC1 and PC2 space
v=zeros(1,length(bfAll.frameTime));
v(:,offset+1:offset+length(behaviorZ))=behaviorZ;
behaviorZ=v;
%hiResBehaviorZ=behaviorZ(diff(BF2stackIdx)>0);
v=zeros(length(CLV),2);
v(offset+1:offset+length(projections),:)=projections;
projections=v;
behaviorTheta=atan2(projections(:,1),projections(:,2));
behaviorTheta=unwrap(behaviorTheta);
behaviorTheta=imfilter(behaviorTheta,filterKernal);
behaviorZTrack=behaviorZ(bfRange);
projectionsTrack=projections(bfRange,:);
behaviorThetaTrack=behaviorTheta(bfRange);



%% make some sort of ethogram;
fcolor=[0 1 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
turncolor=[0 0 1];%[117 112 179]/256;
pausecolor=[255 217 50]/256;
ethocolormap=[bcolor;pausecolor;fcolor;turncolor];
ethoZ=hiResBehaviorZ(hiResRange);
ethoZ=ethoZ-nanmean(ethoZ);
ethoZ(isnan(ethoZ))=0;
Vcluster=smooth(hiResV(hiResRange),100);
idx=sign(Vcluster);
idx(abs(Vcluster)<.00005)=0;
ethogram=idx;
pausing=find(ethogram==0);

idxgmm=kmeans(ethoZ,2,'start',[ 0 5]');
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
ethogram(hiResV(hiResRange)==0)=nan;
ethoTrack=ethogram;
ethoTrack=interp1(hiResFrameTime,ethoTrack,frameTimeTrack,'nearest');

save([dataFolder filesep 'positionData'], 'xPos','yPos','hiResV','behaviorZTrack',...
    'projectionsTrack')
save([dataFolder filesep 'heatData'],'G2','R2','gRaw','rRaw','Ratio2',...
    'gPhotoCorr','rPhotoCorr','acorr','cgIdx','cgIdxRev','ethoTrack','hasPointsTime');
