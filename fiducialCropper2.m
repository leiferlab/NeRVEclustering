%%%%%%%%%%%
% Fiducial cropper 2 takes straightened coordinates form a pointStats file
% and crops out regions corresponding to the unstraighted form. It then
% averages over those ROIs to get the intensity at that time. 
function fiducialCropper2( dataFolder)
if nargin==0
display('Select PointStatsFile');
pointStatsFile=uipickfiles();
pointStatsFile=pointStatsFile{1};
dataFolder=fileparts(pointStatsFile);
display('Select Image Folder');
imageFolder=uipickfiles('filterspec',dataFolder);
imageFolder=imageFolder{1};
else
    pointStatsFile=[dataFolder filesep 'PointStatsNew.mat'];
    psFolder=dir([dataFolder filesep 'CLstraight*']);
    imageFolder=[dataFolder filesep psFolder.name];
end
load(pointStatsFile)

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
    
    rawPointsId=pointStats2(presentIdx(i)).trackIdx(~isnan(pointStats2(presentIdx(i)).trackIdx));
    rawPointsId=(1:length(rawPointsId))';
    

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

[bfTime,hiResFrameTime,hasPoints,bfRange,hiResRange,hasPointsTime,lookup]=dataTimeAlignment(dataFolder,fiducialPoints);

