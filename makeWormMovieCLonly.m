dataFolder='V:\20150617\BrainScanner20150617_143621 - GFP - worm1';
%%
imSize=[1200 600];
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.1 0.01], [0.1 0.01]);

%% load Fiducials file
fiducialFile=dir([dataFolder filesep '*iducial*']);
fiducialFile={fiducialFile.name}';
if length(fiducialFile)~=1
        display('Select model file');

    fiducialFile=uipickfiles('FilterSpec',dataFolder);
    fiducialFile=load(fiducialFile{1});
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;
else
    fiducialFile=load([dataFolder filesep fiducialFile{1}]);
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;

end

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));

%%
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');

%%
behavior=load([dataFolder filesep 'ManualBehavior']);
behavior=behavior.behavior;

hiResBehavior=behavior(diff(BF2stackIdx)>0);


%% load eigenworm data and centerlines
temp=load([dataFolder filesep 'Behavior Analysis' filesep 'eigproj_20141212']);
temp=temp.eigproj_20141212;
offset=9700;


%% load centerline data
centerline=load([dataFolder filesep 'Behavior Analysis' filesep 'centerline']);

if isfield(centerline,'offset')
offset=centerline.offset;
else
    offset=0;
end
wormcentered=centerline.wormcentered;
centerline=centerline.centerline;

%[dangle, f]=centerline2AngleMap(wormcentered);

%%
close all
wormcentered2=smooth2a(wormcentered,5,31);
[wormcenteredx,wormcenteredy]=gradient(wormcentered2);
subplot(3,1,1);imagesc(wormcentered);
subplot(3,1,2);imagesc(sign(wormcenteredy./wormcenteredx))
CLV=sum((wormcenteredy.*wormcenteredx));
subplot(3,1,3);plot(CLV);
xlim([0 length(wormcentered)]);

%%
if offset>0
v= [zeros(1,offset); CLV'];

else
    v=CLV';
end
zeroPad=zeros(size(bfAll.frameTime,1)-length(v),1);
v=[v;zeroPad];
v(1)=0;
CLV=v;
CLbehavior=sign(smooth(-CLV,200));

hiResCLbehavior=interp1(bfAll.frameTime,CLbehavior,hiResData.frameTime);
%%
centerLinePosition=squeeze(mean(centerline,1));
centerLinePosition=bsxfun(@minus,centerLinePosition,mean(centerLinePosition,2));
v=zeros(2,length(behavior));
v(:,offset+1:offset+length(centerLinePosition))=centerLinePosition;
centerLinePosition=v';

hiResCLposition=[interp1(1:length(centerLinePosition),centerLinePosition(:,1),bfIdxLookup),...
    interp1(1:length(centerLinePosition),centerLinePosition(:,2),bfIdxLookup)];

hiResCLposition=hiResCLposition*1/557; % 1mm per 557 pixels
stageCamAngle=81.5;
stageCamAngle=stageCamAngle*pi/180;
%rotation matrix that takes motion of stage direction to motion in low mag
%behavior image
Rmatrix=[-cos(stageCamAngle) sin(stageCamAngle);...
    -sin(stageCamAngle) -cos(stageCamAngle)];
hiResCLposition=(Rmatrix'*hiResCLposition')';


%% load eigen behavior
temp=load([dataFolder filesep 'Behavior Analysis' filesep 'eigproj_20141212']);
temp=temp.eigproj_20141212;

hasPoints=cellfun(@(x) ~isempty(x{1}), fiducialPoints,'uniformoutput',0);
hasPoints=find(cell2mat(hasPoints));
nTimes=length(hasPoints);


hasPointsTime=hiResData.frameTime(diff(hiResData.stackIdx)==1);
hasPointsTime=hasPointsTime(hasPoints);
hasPointsTime=hasPointsTime-min(hasPointsTime);
frameTime=hiResData.frameTime;

hiResRangeTrack=(find(diff(hiResData.stackIdx)>0));
hiResRangeTrack=hiResRangeTrack(hasPoints);
hiResRange=hiResRangeTrack(1):hiResRangeTrack(end);
hiResFrameTime=frameTime(hiResRange);
timeStart=min(hiResFrameTime);
hiResFrameTime=hiResFrameTime-timeStart;

bfRange=find(diff(BF2stackIdx)==1);
bfRange=bfRange(hasPoints(1)):bfRange(hasPoints(end));
bfTime=bfAll.frameTime(bfRange);
bfTime=bfTime-bfTime(1);

eigenBehavior=zeros(size(temp,1),length(behavior));
eigenBehavior(:,offset+1:offset+length(temp))=temp;

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

v=zeros(1,length(behavior));
v(:,offset+1:offset+length(behaviorZ))=behaviorZ;
behaviorZ=v;
hiResBehaviorZ=behaviorZ(diff(BF2stackIdx)>0);
v=zeros(length(behavior),2);
v(offset+1:offset+length(projections),:)=projections;
projections=v;
behaviorTheta=atan2(projections(:,1),projections(:,2));
hiResBehaviorTheta=behaviorTheta(diff(BF2stackIdx)>0);
hiResEigen=projections(diff(BF2stackIdx)>0,:);

%hiResV=gradient(D)./gradient(hiResData.frameTime);

%%


xPos=hiResData.xPos/10000;
%xPos([0; diff(xPos)]==0)=nan;
xPos=inpaint_nans(xPos);
yPos=hiResData.yPos/10000;
%yPos([0; diff(yPos)]==0)=nan;
yPos=inpaint_nans(yPos);
xPos=xPos-1*hiResCLposition(:,2);
yPos=yPos-1*hiResCLposition(:,1);


xPos=smooth(xPos,100);
yPos=smooth(yPos,100);

hiResxPos=xPos(hiResRange);
hiResyPos=yPos(hiResRange);
hiResxPos=hiResxPos-(max(hiResxPos)+min(hiResxPos))/2;
hiResyPos=hiResyPos-(max(hiResyPos)+min(hiResyPos))/2;


figure
h=colorplot(hiResxPos(1:100:end),hiResyPos(1:100:end),hiResFrameTime(1:100:end));
axescenter
colorbar
set(gca,'fontsize',15);

%%
frameTimeTrack=hiResData.frameTime((find(diff(hiResData.stackIdx)>0)));
frameTimeTrack=frameTimeTrack(hasPoints);
frameTimeTrack=frameTimeTrack-min(frameTimeTrack);

xPosTrack=hiResData.xPos((find(diff(hiResData.stackIdx)>0)));
xPosTrack=xPosTrack(hasPoints);
yPosTrack=hiResData.yPos((find(diff(hiResData.stackIdx)>0)));
yPosTrack=yPosTrack(hasPoints);




v=[diff(xPos) diff(yPos)];
v=sqrt(sum(v.^2,2));
D=[0;cumsum(v)];
D=smooth(D,100);

hiResV=gradient(D)./gradient(hiResData.frameTime);
hiResV=hiResV.*hiResCLbehavior;

%% make ethogram
fcolor=[0 1 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
pausecolor=[0 0 1];%[117 112 179]/256;

idx=kmeans(smooth(hiResV(hiResRange),100),3,'start',[-.2 0 .2]');
idx=idx-2;
ethogram=idx;
pausing=ethogram==0;

 pausing=abs(smooth(pausing,100))>.94;

ethogram(ethogram==0 & ~pausing)=nan;
ethogram=colNanFill(ethogram);

ethogram(hiResV(hiResRange)==0)=nan;
ethoTrack=ethogram;
ethoTrack(hiResV(hiResRange)==0)=nan;
ethoTrack=interp1(hiResFrameTime,ethoTrack,frameTimeTrack);


%%
hiResActivityFolder=[dataFolder filesep 'fiducials' filesep 'hiResActivityFolder3Dtest'];
hiResSegmentFolder=[dataFolder filesep 'fiducials' filesep  'hiResSegmentFolder3Dtest'];
lookupFolderX=[dataFolder filesep 'fiducials' filesep  'lookupFolderXtest'];
lookupFolderY=[dataFolder filesep 'fiducials' filesep  'lookupFolderYtest'];
metaFolder=[dataFolder filesep 'fiducials' filesep  'metaDataFolder3Dtest'];
hiResSegmentFolderRaw=[dataFolder filesep 'fiducials' filesep  'hiResSegmentFolder3Dtest_raw'];
hiResActivityFolderRaw=[dataFolder filesep 'fiducials' filesep  'hiResActivityFolder3Dtest_raw'];
hiResSegmentFolderS1=[dataFolder filesep 'fiducials' filesep  'hiResSegmentFolder3Dtest_S1'];
hiResActivityFolderS1=[dataFolder filesep 'fiducials' filesep  'hiResActivityFolder3Dtest_S1'];
lowResFolder=[dataFolder filesep 'fiducials' filesep  'lowResFluor'];
lowResBFFolder=[dataFolder filesep 'fiducials' filesep  'lowResBF'];
%imageFolder2=[dataFolder filesep 'imageMovieAllBalls'];
imageFolder=[dataFolder filesep 'rawVideo2'];

mkdir(imageFolder);
%mkdir(imageFolder2);
%%
stackList=1:973;
[XMAX] = [max(xPos)];
XMIN=[min(xPos)];
[YMAX] = max(yPos);
[YMIN] = min(yPos);
colorLookup=linspace(min(behaviorZ),max(behaviorZ),256);
colorMap=jet(256);


if (XMAX-XMIN)>(YMIN-YMAX)
    range=XMAX-XMIN;
    ycenter=(YMIN+YMAX)/2;
    yrange=ycenter+[-range/2, range/2];
    xrange=[XMIN XMAX];
else
        range=YMAX-YMIN;
    xcenter=(XMIN+XMAX)/2;
    xrange=xcenter+[-range/2, range/2];
    yrange=[YMIN YMAX];
    
end


%%
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
S2AHiRes=alignments.S2AHiRes;
lowResFluor2BF=alignments.lowResFluor2BF;

    rect1=S2AHiRes.rect1;
    rect2=S2AHiRes.rect2;
    t_concord=S2AHiRes.t_concord;
    Rsegment=S2AHiRes.Rsegment;
    padRegion=S2AHiRes.padRegion;
    
   lowResRsegment=lowResFluor2BF.Rsegment;
   lowResT=lowResFluor2BF.t_concord;

%%
fclose all;

    Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat'] );
frameRange=find(hiResData.stackIdx==stackList(1),1,'first'):find(hiResData.stackIdx==stackList(end),1,'last');
frameRange=68000:100000;
%% prep vidobj for avi files and .dat file

%search in dataFolder for avi's without the HUDS, one with the fluor and
%one without.
aviFiles=dir([dataFolder filesep '*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
if length(aviFiles)==2
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviorMovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluorMovie=[dataFolder filesep aviFiles{aviFluorIdx}];
else
    display('Select avi files, behavior and then low mag fluor');
    movies=uipickfiles('FilterSpec',dataFolder);
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);


%%
pointHistory=[];
frameCounter=1;
imageFolder=[dataFolder filesep 'rawVideoFeeds2'];
mkdir(imageFolder);
startTime=hiResData.frameTime(frameRange(1));

for iStack=1:10:length(frameRange)
    tic
  
        %%
        try
    hiResIdx=frameRange(iStack);
    %imName=['image' num2str(stackIdx,'%3.5d') '.tif'];
    %matName=['image' num2str(iStack,'%3.5d') '.mat'];
%     metaData=load([metaFolder filesep matName]);
%     metaData=metaData.metaData.metaData;
 %   wormInfo=imfinfo([hiResSegmentFolder filesep imName]);
%     
% worm=stackLoad([hiResSegmentFolder filesep imName]);
% activity=stackLoad([hiResActivityFolder filesep  imName]);
% if metaData.zVoltage(1)>metaData.zVoltage(end);
%     worm=flipdim(worm,3);
%     activity=flipdim(activity,3);
% end
% worm=normalizeRange(worm);
% activity=normalizeRange(activity);
%%

bfFrame=round(bfIdxLookup(hiResIdx));
bfImage=read(behaviorVidObj,bfFrame);
bfImage=squeeze(bfImage(:,:,1,:));
bfImage=pedistalSubtract(double(bfImage));

bfCenterline=centerline(:,:,bfFrame-offset);

%c
    filename=[imageFolder filesep 'frame' num2str(iStack,'%3.5d')];
 %   filename2=[imageFolder2 filesep 'frame' num2str(frameCounter,'%3.5d')];
   % worm
     imagesc(bfImage);
        text(150, 150, {'Behavior (10x)'},'fontsize',12,'color',[1 1 1]);

    colormap gray
    axis equal off tight
    hold on
   plot( bfCenterline(:,2),bfCenterline(:,1));
hold off
text(150,950,['t = ' num2str(hiResData.frameTime(hiResIdx)-startTime,'%3.1f') 's'],'color','w','fontsize',20)

%eigenPosture
  
set(gcf,'color',[1 1 1])
drawnow
saveas(gcf,filename,'jpeg');


    display(['completed stack:' num2str(iStack) ' in ' num2str(toc) 's']);
        catch ME
      display(['error stack:' num2str(iStack) ' in ' num2str(toc) 's']);
ME
      %save(['errorStack' num2str(iStack)], 'ME') 
        end
        
end
