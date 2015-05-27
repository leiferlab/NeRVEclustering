dataFolder='F:\20141212\BrainScanner20141212_145951\';
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
centerline=load([dataFolder filesep 'Behavior Analysis' filesep 'centerline_4245']);
centerline=centerline.centerline_4245;

[dangle, f]=centerline2AngleMap(centerline);
v=zeros(1,length(behavior));
v(:,offset+1:offset+length(f.v))=f.v;
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
xPos([0; diff(xPos)]==0)=nan;
xPos=inpaint_nans(xPos);
yPos=hiResData.yPos/10000;
yPos([0; diff(yPos)]==0)=nan;
yPos=inpaint_nans(yPos);
xPos=xPos-1*hiResCLposition(:,2);
yPos=yPos-1*hiResCLposition(:,1);
filterKernal=gausswin(300);
filterKernal=filterKernal/sum(filterKernal);
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
fcolor=[0 .7 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
turncolor=[0 0 1];%[117 112 179]/256;
pausecolor=[.2 .2 .2];

load([dataFolder filesep 'ethoTrack'])

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
imageFolder=[dataFolder filesep 'imageMovieAllBalls4'];
%imageFolder2=[dataFolder filesep 'imageMovieAllBalls'];

mkdir(imageFolder);
%mkdir(imageFolder2);
%%
stackList=500:973;
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
masterData=load([dataFolder filesep 'refStackStraight']);
fiducials=masterData.masterFiducials;
data=load([dataFolder  filesep 'heatData20141212']);

%%
figure
possibleB=[2    61    67];
possibleF=    [20    18    60    10] ;%old result:[20    18    42    10];
possibleP=[14    19    35    29    33    45     7];%[  19    35    29    64    33    45     7];

colorBalls=.3+zeros(size(fiducials));
colorBalls(possibleP,3)=1;
colorBalls(possibleB,1)=1;

transp=.2*ones(1,length(fiducials));
transp([possibleB possibleF possibleP])=.8;


colorBalls(possibleF,2)=1;
scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',.002,'color',colorBalls,'transp',transp);axis equal off tight

G2=data.G2;
G2=colNanFill(G2')';

G2ColorsLookup=G2;
G2ColorsLookup=G2ColorsLookup/1.25;
G2ColorsLookup(G2ColorsLookup>1)=1;
G2ColorsLookup=colNanFill(G2ColorsLookup')';
G2ColorsLookup=ceil(G2ColorsLookup*100);

G2ColorsLookup(isnan(G2ColorsLookup))=1;
G2ColorsLookup(G2ColorsLookup<1)=1;
Pmap=cbrewer('seq','Blues',100);
Bmap=cbrewer('seq','YlOrRd',100);
Fmap=cbrewer('seq','YlGnBu',100);
Fmap=parula(100);
%%
pointHistory=[];
frameCounter=1;

for iStack=1:length(stackList)
    tic
  
        %%
        try
    stackIdx=stackList(iStack);
    imName=['image' num2str(stackIdx,'%3.5d') '.tif'];
    matName=['image' num2str(iStack,'%3.5d') '.mat'];
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
frameRange=find(hiResData.stackIdx==stackIdx);
bfFrame=round(bfIdxLookup(frameRange));
[bfFrameUnique,ia]=unique(bfFrame);
bfImage=read(behaviorVidObj,round(bfFrameUnique([1,end])));
bfImage=squeeze(bfImage(:,:,1,:));

bfeigenProj=projections(bfFrameUnique,:);
bfeigenZ=behaviorZ(bfFrameUnique);
bfCenterline=centerline(:,:,bfFrameUnique-offset);
hiResApproxIdx=frameRange;

xPosFrame=interp1(xPos,hiResApproxIdx);
yPosFrame=interp1(yPos,hiResApproxIdx);

%bfIdx=round(linspace(1,size(bfImage,3),length(wormInfo)));
bfIdx=1:size(bfImage,3);

behaviorState=ethoTrack(iStack);

switch behaviorState
    case -1
        color=bcolor;
        stateText='Back';
    case 4
        color = turncolor;
        stateText='Turn';

    case 1
        color = fcolor;
        stateText='Forward';
    case 0
        color = pausecolor;
        stateText='Pause';        

end
%%
  colorT=G2ColorsLookup(:,iStack);

%    colorBalls=.3+zeros(size(fiducials));
% colorBalls(possibleP,:)=Fmap(G2ColorsLookup(possibleP,iStack),:);
% colorBalls(possibleB,:)=Fmap(G2ColorsLookup(possibleB,iStack),:);
% colorBalls(possibleF,:)=Fmap(G2ColorsLookup(possibleF,iStack),:);
colorBalls=Fmap(G2ColorsLookup(:,iStack),:);
ballSize=G2(:,iStack);
%ballSize=ballSize-nanmin(G2(:));
ballSize(ballSize<0)=0;
ballSize(isnan(ballSize))=0;
ballSize=ballSize*.0015+.5*.002;
transp=.6*ones(1,length(fiducials));
%transp([possibleB possibleF possibleP])=.8;
subplot(2,2,[1 2]);
cla

scatter3sph(fiducials(:,1)/4000,fiducials(:,2)/4000,3*fiducials(:,3)/4000,...
    'size',ballSize,'color',colorBalls,'transp',transp);axis equal off tight

colormap(Fmap)
%h=lcolorbar(0:.5:2);
caxis([0 1.3])
%set(h,'Location','westOutside')
nBalls=5;
midBall=round(nBalls/2);
refY=linspace(.0825,.105, nBalls);
refX=ones(1,nBalls)*.155;
refZ=zeros(1,nBalls);
refI=linspace(0,1.333,nBalls);
refSize=refI*.0015+.5*.002;

refColor=interp1(Fmap,linspace(1,100,nBalls));


hold on

scatter3sph(refX,refY,refZ,...
    'size',refSize,'color',refColor);axis equal off tight




sph=gca;
 htext=text(0.025,0.12,0,stateText,'parent',sph,'Color',color,'FontSize',20);
 htext2=text(refX+.005,refY,refZ,strread(num2str(refI,'%6.2g'),'%s'),'parent',sph,'Color',[0 0 0],'FontSize',12);
 htext3=text(refX(midBall)+.02,refY(midBall),refZ(midBall),'\DeltaF/F0'...
     ,'parent',sph,'Color',[0 0 0],'FontSize',16,'Rotation',-90,'HorizontalAlignment','center');
% ylimits=[.079, .113];
% xlimits=[.03,.18];
%  xlim(xlimits);
% ylim(ylimits);


 % h=colorbar('peer',sph);
% 
% h.ColormapMode='manual';
% htext=text(0.05,0.11,0,stateText,'parent',sph,'Color',color,'FontSize',20);
% 
% set(get(h,'YLabel'),'String','\Delta F/F0')
view([0,90])

freezeColors
%  if iStack==1
%    %  hrect=getrect();
%      hrect(2)=max(hrect(2),ylimits(1));
%     hrect(4)=min(hrect(4)+hrect(2),ylimits(2))-hrect(2);
%      hrect(1)=max(hrect(1),xlimits(1));
%      hrect(3)=min(hrect(3),xlimits(2)-hrect(1));
%           [hdrawx,hdrawy]=getpts();
%           hdrawx=smooth(hdrawx,3);
%           hdrawy=smooth(hdrawy,3);
%      ylimits(1)=min(hdrawy);
%      ylimits(2)=max(hdrawy);
%      xlimits(1)=min(hdrawx);
%      
% 
%  end 
 

 xlimits(2)=.18;
 xlim(xlimits);
ylim(ylimits);
 rectangle('Position',hrect)
plot(hdrawx,hdrawy,'black','linewidth',2)
%%
for iSlice=1:length(bfIdx);
    filename=[imageFolder filesep 'frame' num2str(frameCounter,'%3.5d')];
%    filename2=[imageFolder2 filesep 'frame' num2str(frameCounter,'%3.5d')];
    %worm
%     subplot(2,2,1)
%     
%     imagescwithnan(worm(:,:,iSlice),...
%     hot,[0 0 0 ]);
% colorbar off
%    % imagesc(worm(:,:,iSlice));
%     set(gcf,'color','w');
%     caxis([0 1]);
%     freezeColors
% axis equal off tight
% caxis([0,1]);
%behavior
subplot(2,2,3)
    imagesc(bfImage(:,:,bfIdx(iSlice)));
    colormap gray
    freezeColors
    axis equal off tight
    hold on
   plot( bfCenterline(:,2,bfIdx(iSlice)),bfCenterline(:,1,bfIdx(iSlice)));

   text(136,136,'Behavior','Color',[1 1 1],'FontSize',12)
   hold off
%eigenPosture
%     subplot(2,2,1);
% %delete(gca)
% colorIdx=interp1(colorLookup,1:256,bfeigenZ(bfIdx(iSlice)),'nearest');
%     colorVector=colorMap(colorIdx,:);
%     scatter3(bfeigenProj(bfIdx(iSlice),1),bfeigenProj(bfIdx(iSlice),2),...
%         bfeigenZ(bfIdx(iSlice)),'markerFacecolor',color,'markeredgecolor','none');        
%    xlim([-8,8]);
%     ylim([-8,8]);
%     zlim([-20 20]);
%     xlabel('Eig1');ylabel('Eig2');zlabel('Eig3');
%     
%        imagescwithnan(activity(:,:,iSlice),...
%     pmkmp(123,'LinearL'),[0 0 0 ]);
% colorbar off
%    % imagesc(worm(:,:,iSlice));
%     set(gcf,'color','w');
%     caxis([0 1]);
%     freezeColors
% axis equal off tight
% caxis([0,1]);
% 

    
    
    %position in space
    subplot(2,2,4)
  
    scatter(xPosFrame(bfIdx(iSlice)),yPosFrame(bfIdx(iSlice)),'r',...
        'markerFacecolor',color,'markeredgecolor','none');        


    axis equal 
    
    text(xrange(1)+diff(xrange)/8,yrange(2)-diff(yrange)/8,'Position','fontsize',12);
        xlim(xrange);
    ylim(yrange);
    xlabel('x (mm)');
%     hold on
% pointHistory=[pointHistory;[xPosFrame(bfIdx(iSlice)),xPosFrame(bfIdx(iSlice))]];
%  historyLength=size(pointHistory,1);
% plot(pointHistory(1:10:end,1)...
%     ,pointHistory(1:10:end,2))
% hold off


set(gcf,'color',[1 1 1])
drawnow
saveas(gcf,filename,'jpeg');
frameCounter=frameCounter+1;
end


    display(['completed stack:' num2str(iStack) ' in ' num2str(toc) 's']);
        catch ME
      display(['error stack:' num2str(iStack) ' in ' num2str(toc) 's']);
save(['errorStack' num2str(iStack)], 'ME') 
        end
        
end
