%dataFolder=uipickfiles;
%dataFolder=dataFolder{1};
dataFolder=uipickfiles;
dataFolder=dataFolder{1};
%%
%'F:\20141212\BrainScanner20141212_145951\';
imSize=[1200 600];
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
%% load data
wormFiducialIntensityFile=dir([dataFolder filesep 'wormFiducialIntensities*.mat']);
load([dataFolder filesep wormFiducialIntensityFile.name])

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear','extrap');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime),'linear','extrap');

%% load Fiducials file
 [fiducialPoints,z2ImageIdxOffset]=loadFiducialPoints(dataFolder);
%%
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest','extrap');





%% load manual behavior
% behavior=load([dataFolder filesep 'ManualBehavior']);
% behavior=behavior.behavior;
% 
% hiResBehavior=interp1(1:length(behavior),behavior,bfIdxLookup);
% hiResBehavior(abs(hiResBehavior)~=1)=nan;
% hiResBehavior=inpaint_nans(hiResBehavior);
% hiResBehavior=sign(hiResBehavior);
% offset=9700;


%%
hasPoints=cellfun(@(x) size(cell2mat(x),1), fiducialPoints,'uniformoutput',0);
hasPoints=cell2mat(hasPoints);
hasPoints=find(hasPoints>max(hasPoints)*.9);
hasPoints=hasPoints(hasPoints<max(BF2stackIdx));

%%
bfRange=[1 find(diff(BF2stackIdx)==1)];
hasPoints=hasPoints(hasPoints<length(bfRange));
hasPoints=min(hasPoints):max(hasPoints);

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


bfRange=bfRange(hasPoints(1)):bfRange(hasPoints(end));
bfTime=bfAll.frameTime(bfRange);
bfTime=bfTime-bfTime(1);


%% load centerline data
centerlineFile=load([dataFolder filesep 'BehaviorAnalysis' filesep 'centerlines']);
fieldNames=fieldnames(centerlineFile);
eigIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'eig'));
for iField=1:length(fieldNames)
    currentStr=fieldNames{iField};
if strfind(currentStr,'wormcen')
    wormcentered=centerlineFile.(currentStr);
elseif strfind(currentStr,'line')
    centerline=centerlineFile.(currentStr);
end
end

if isfield(centerlineFile,'offset');
offset=centerlineFile.offset;
else
    offset=0;
end

if any(eigIdx)
    eigenProj=centerlineFile.(fieldNames{eigIdx});
else
    eigenProj=load([dataFolder filesep 'BehaviorAnalysis' filesep 'eigenProj']);
    fieldNames=fieldnames(eigenProj);
eigenProj=eigenProj.(fieldNames{1});

    
end

    
    

    


%[dangle, f]=centerline2AngleMap(wormcentered);

%%
close all

wormcentered2=smooth2a(wormcentered,5,31);
[wormcenteredx,wormcenteredy]=gradient(wormcentered2);
CLV=sum((wormcenteredy.*wormcenteredx));

subplot(3,1,1);imagesc(wormcentered);
subplot(3,1,2);imagesc(sign(wormcenteredy./wormcenteredx))
subplot(3,1,3);plot(CLV);
xlim([0 length(wormcentered)]);

%%
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

%%
centerLinePosition=squeeze(mean(centerline,1));
centerLinePosition(centerLinePosition==0)=nan;
centerLinePosition=inpaint_nans(centerLinePosition);
centerLinePosition=bsxfun(@minus,centerLinePosition,mean(centerLinePosition,2));
v=zeros(2,length(CLV));
v(:,offset+1:offset+length(centerLinePosition))=centerLinePosition;
zeroPad=zeros(size(bfAll.frameTime,1)-length(v),1);
v=[v;zeroPad];

centerLinePosition=v';



hiResCLposition=[interp1(1:length(centerLinePosition),centerLinePosition(:,1),bfIdxLookup,'nearest','extrap'),...
    interp1(1:length(centerLinePosition),centerLinePosition(:,2),bfIdxLookup,'nearest','extrap')];

hiResCLposition=hiResCLposition*1/557; % 1mm per 557 pixels
stageCamAngle=90;
stageCamAngle=stageCamAngle*pi/180;
%rotation matrix that takes motion of stage direction to motion in low mag
%behavior image
Rmatrix=[-cos(stageCamAngle) sin(stageCamAngle);...
    -sin(stageCamAngle) -cos(stageCamAngle)];
hiResCLposition=(Rmatrix'*hiResCLposition')';
filterKernal=gausswin(50);
filterKernal=filterKernal-min(filterKernal(:));
filterKernal=filterKernal/sum(filterKernal);

%% load eigen behavior


% bfRange=find(diff(BF2stackIdx)>0);
% bfRange=bfRange(hasPoints(1)):bfRange(hasPoints(end));
% 
% bfTime=bfAll.frameTime(bfRange);
% bfTime=bfTime-bfTime(1);

eigenBehavior=zeros(size(eigenProj,1),length(CLV));
eigenBehavior(:,offset+1:offset+length(eigenProj))=eigenProj;
%v=zeros(size(eigenBehavior,1),length(bfAll.frameTime));
%v(:,offset+1:offset+length(eigenBehavior))=eigenBehavior;
%zeroPad=zeros(size(bfAll.frameTime,1)-length(v),1);
%v=[v;zeroPad];
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

% %%
% projections=eigenBehavior(1:2,:)';
% behaviorZ=eigenBehavior(3,:)';

%%
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
hiResBehaviorTheta=behaviorTheta(diff(BF2stackIdx)>0);
hiResEigen=projections(diff(BF2stackIdx)>0,:);






%%
subplot(2,1,1);
plotyy(1:length(xPosStage),xPosStage,1:length(xPosStage),yPosStage);
subplot(2,1,2);
plot(hiResCLposition);
%%
xPosStage=hiResData.xPos/10000;
xPosStage([0; diff(xPosStage)]==0)=nan;
xPosStage=inpaint_nans(xPosStage);
yPosStage=hiResData.yPos/10000;
yPosStage([0; diff(yPosStage)]==0)=nan;
yPosStage=inpaint_nans(yPosStage);
xPos=xPosStage+1*hiResCLposition(:,2);
yPos=yPosStage+1*hiResCLposition(:,1);
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
hiResxPos=hiResxPos-(max(hiResxPos)+min(hiResxPos))/2;
hiResyPos=hiResyPos-(max(hiResyPos)+min(hiResyPos))/2;


figure
%h=colorplot(hiResxPos(1:100:end),hiResyPos(1:100:end),hiResFrameTime(1:100:end));

scatter(hiResxPos,hiResyPos,[],hiResFrameTime,'.');
scalebar= [max(hiResxPos) min(hiResyPos)];
scalebar=[scalebar ; scalebar+[1 0]];
hold on
plot(scalebar(:,1)-1,scalebar(:,2),'black','LineWidth',2)
save([dataFolder filesep 'positionData'],'hiResxPos','hiResyPos','hiResFrameTime');
colormap parula
axis equal
%axescenter
%colorbar
axis off
set(gca,'fontsize',15);

%%

filterKernal2=gausswin(100);
filterKernal2=filterKernal2-min(filterKernal2(:));
filterKernal2=filterKernal2/sum(filterKernal2);
filterFactor2=imfilter(ones(size(xPos)),filterKernal2);


frameTimeTrack=hiResData.frameTime((find(diff(hiResData.stackIdx)>0)));
frameTimeTrack=frameTimeTrack(hasPoints);
frameTimeTrack=frameTimeTrack-min(frameTimeTrack);

xPosTrack=xPos((find(diff(hiResData.stackIdx)>0)));
xPosTrack=xPosTrack(hasPoints);
yPosTrack=yPos((find(diff(hiResData.stackIdx)>0)));
yPosTrack=yPosTrack(hasPoints);


 xV=imfilter((xPos),filterKernal2)./filterFactor2;
 yV=imfilter((yPos),filterKernal2)./filterFactor2;


v=[gradient(xV) gradient(yV)];
v=sqrt(sum(v.^2,2));
D=[cumsum(v)];

hiResV=gradient(D)./median(gradient(hiResData.frameTime));
hiResV=hiResV.*hiResCLbehavior;
hiResV=smooth(hiResV,500);

figure;
plot(hiResData.frameTime(hiResRange)-hiResData.frameTime(hiResRange(1)),hiResV(hiResRange),...
    'color',[106 61 154]/256,'linew',2);
axis tight;
ylim([-.15, .15])
xlabel('Time (s)');
ylabel('Velocity (mm/s)');
set(gca,'box','off');
%xaxescenter
set(gca,'fontSize',15);



%%
%reversalTimes=[5:13 41:57 81:87 110:120 160:170];
%behaviorTrack=hiResStackBehavior(hasPoints);
behaviorZtrack=hiResBehaviorZ(hasPoints);
behaviorThetaTrack=hiResBehaviorTheta(hasPoints);
projectionsTrack=hiResEigen(hasPoints,:);
%vTrack=hiResV(hasPoints);


%%
figure;
h(1)=subplot(3,1,1);
lineColor=[.3 .6 .3];
%colorplot(bfTime,projections(bfRange,1),bfTime);
plot(bfTime,projections(bfRange,1),'color',lineColor,'linew',2);
axis tight
ylim([-8,8])

xaxescenter
ylabel('Eig1')
set(gca,'fontSize',15);

h(2)=subplot(3,1,2);
%colorplot(bfTime,projections(bfRange,2),bfTime);
plot(bfTime,projections(bfRange,2),'color',lineColor,'linew',2);


ylim([-8,8])
xaxescenter
axis tight
ylabel('Eig2')
set(gca,'fontSize',15);

h(3)=subplot(3,1,3);
%colorplot(bfTime,behaviorZ(bfRange),bfTime);
plot(bfTime,behaviorZ(bfRange),'color',lineColor,'linew',2);

ylim([-8,8])
axis tight
ylabel('Eig3');
xlabel('Time(s)');
ylim([-30,30])
xaxescenter
set(gca,'fontSize',15);
linkaxes(h,'x')

colormap parula


%% make some sort of ethogram;
fcolor=[0 1 0];%[27 158 119]/256;
bcolor=[1 0 0];%[217 95 2]/256;
turncolor=[0 0 1];%[117 112 179]/256;
pausecolor=[255 217 50]/256;
ethocolormap=[bcolor;pausecolor;fcolor;turncolor];

% ethoV=interp1(hiResFrameTime,hiResV(hiResRange),1:.01:max(hiResFrameTime));
% ethoB=interp1(hiResFrameTime,hiResBehavior(hiResRange),1:.01:max(hiResFrameTime));
load('Y:\CommunalCode\3dbrain\Behavior\eigenProj3GMM','gmm')


% obj=gmdistribution.fit((hiResV(hiResRange)),3,'start',...
%     2+round(smooth(hiResBehavior(hiResRange),100)));
% temp=hiResV;
% temp(abs(temp)<.02)=0;
% obj=gmdistribution.fit((hiResV(hiResRange)),3,'start',...
%     2+sign(temp(hiResRange)));
% 
% [idx,nlogl,P,logpdf,M] = cluster(obj,smooth(hiResV(hiResRange),100));
% 
% [~,ia]=sort(obj.mu,'ascend');
%ethoZ=interp1(bfTime,behaviorZ(bfRange),(hiResFrameTime));
ethoZ=hiResBehaviorZ(hiResRange);
ethoZ=ethoZ-mean(ethoZ);
Vcluster=smooth(hiResV(hiResRange),100);
idx=kmeans(Vcluster,3,'start',[-.2 0 .2]');
idx=idx-2;
idx(Vcluster>.1)=1;
idx(Vcluster<-0)=-1;
idx(idx~=1 & idx~=-1)=0;
ethogram=idx;
pausing=find(ethogram==0);

idxgmm=kmeans(ethoZ,2,'start',[ 0 5]');
%%
ethogram(((abs(ethoZ)>2*std(ethoZ))& (ethogram>=0)))=2;
%ethogram((abs(ethoZ)>2*std(ethoZ)))=2;


%pausing=(ethogram==0);

%% Kill off short behaviors unless they are reversals
for iBehavior=[0 1 2];
cpause=bwconncomp(ethogram==iBehavior);
shortPause=cellfun(@(x) length(x), cpause.PixelIdxList);
 shortPause= shortPause'<500;
 shortPause=cpause.PixelIdxList(shortPause);
 shortPause=cell2mat(shortPause');
 ethogram(shortPause)=nan;
ethogram=colNanFill(ethogram,'nearest');
end
%%

% ethogram=ethogram(hiResRange);
figure
imagesc(hiResFrameTime,1,ethogram');
colormap(ethocolormap)
caxis([-1,2])
%xlim([0 min(max(hiResFrameTime(hiResV(hiResRange)~=0)),max(hasPointsTime(sum(isnan(G2))<50))')]);
set(gca,'YTick',[]);
ethogram(hiResV(hiResRange)==0)=nan;
ethoTrack=ethogram;
ethoTrack(hiResV(hiResRange)==0)=nan;
ethoTrack=interp1(hiResFrameTime,ethoTrack,frameTimeTrack,'nearest');
%save([dataFolder filesep 'ethoTrack'],'ethoTrack');

%%
centerlineShowFrame=hasPoints(1:round(length(hasPoints)/5):end);
centerlineShowFrame=stack2BFidx(centerlineShowFrame)-offset;
centerlineShowFrame=min(round(centerlineShowFrame),size(centerline,3));
subCenterlines=centerline(:,:,centerlineShowFrame);
figure
for i=1:size(subCenterlines,3)
    temp=subCenterlines(:,:,i);
    temp=flipud(temp);
    temp(:,1)=-temp(:,1);
    temp=bsxfun(@minus,temp,(min(temp)+max(temp))/2);
    temp(:,1)=temp(:,1)+(i-1)*500+250;
    temp(:,2)=temp(:,2)+500;
    
    
    plot(temp(:,1),temp(:,2),'black','linew',2);
        hold on

    scatter(temp(1,1),temp(1,2),'r')
end
timeVec=1:(i)*500;
%imagesc(hiResFrameTime,1,ethogram');

% h=colorplot(hiResFrameTime/max(hiResFrameTime)*2500,1,ethogram);
% colormap(ethocolormap)

%h=colorplot(hiResFrameTime/max(hiResFrameTime)*2500,1,hiResFrameTime);

%set(h,'linew',50);

    axis equal tight off
    

%%
ballSearch=7;
cropSize=size(Rcrop);
vols=cropSize(5);
RvalsAll=zeros(vols,cropSize(4));
GvalsAll=RvalsAll;

center=ceil(cropSize/2);
center=center(1:3);
centerIdx=sub2ind(cropSize(1:3),center(1),center(2),center(3));
[ballx,bally]=meshgrid(1:cropSize(1),1:cropSize(2));
ballx=ballx-center(1);
bally=bally-center(2);
ballr=sqrt(ballx.^2+bally.^2);
ballr=repmat(ballr<ballSearch,1,1,cropSize(3));
progressbar(0,0);
for iVol=1:vols
    for i=1:cropSize(4)
      %  iTime=hasPoints(i);
      iTime=i;
        progressbar(iVol/vols,iTime/nTimes);
        subVolumeR=(Rcrop(:,:,:,iTime,iVol))-90;
     %           subVolumeR=pedistalSubtract(Rcrop(:,:,:,iTime,iVol));

        subVolumeG=(Gcrop(:,:,:,iTime,iVol))-90;
        subVolumeG((subVolumeG)<=0)=nan;
        subVolumeR(subVolumeR<0)=nan;

        subBW=(subVolumeR)>(.5*subVolumeR(center(1),center(2),center(3)));
        subBW=subBW.*ballr;
        cc=bwconncomp(subBW,6);
                    nnz(subBW)

        centerPixPresent=cellfun(@(x) any(x==centerIdx),cc.PixelIdxList);
        if any(centerPixPresent)
            cc.PixelIdxList=cc.PixelIdxList(centerPixPresent);
            areaIdx=1;
        else
        
        nPixRegions=cell2mat(cellfun(@(x) length(x), cc.PixelIdxList,'uniformOutput',0));
        [~,areaIdx]=max(nPixRegions);
        end
        Rvals=subVolumeR(cell2mat(cc.PixelIdxList(areaIdx)'));
        Gvals=subVolumeG(cell2mat(cc.PixelIdxList(areaIdx)'));
        
        RvalsAll(iVol,i)=nanmean(Rvals);
        GvalsAll(iVol,i)=nanmean(Gvals);
    end
end

%%
       RvalsAll(RvalsAll==0|RvalsAll>(nanstd(RvalsAll(:))*10))=nan;

GvalsAll(GvalsAll==0|GvalsAll>(nanstd(GvalsAll(:))*10))=nan;
rejects=(mean(isnan(RvalsAll),2)>.1);
rejectCorrection=cumsum(rejects);
rejectCorrection(rejects)=[];
rejectCorrection=rejectCorrection+(1:length(rejectCorrection))';
GvalsAll(rejects,:)=[];

RvalsAll(rejects,:)=[];

%%

%A=(G2smooth./R2smooth)';
%correct photobleaching red
photoBleachingR=zeros(size(RvalsAll));
photoBleachingG=zeros(size(RvalsAll));

%Fexponent =@(x,xdata) x(1)*exp(x(2)*xdata)+x(3);
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});
fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1001,0,10000];

progressbar(0)
for i=1:size(RvalsAll,1)
    
    try
        %%
     
    progressbar(i/size(RvalsAll,1));
    xVals=(1:size(RvalsAll,2))';
    
    present=(~isnan(RvalsAll(i,:)+GvalsAll(i,:))') ;
    xVals=xVals(present);
    rVals=RvalsAll(i,:)';
    minWindow=150;
    gVals=GvalsAll(i,:)';
        gVals=gVals(present);
    rVals=rVals(present);
    gVals=ordfilt2(gVals,30,true(minWindow,1));
    rVals=ordfilt2(rVals,30,true(minWindow,1));
    gVals(1:(minWindow/2))=gVals((minWindow/2));
    rVals(1:(minWindow/2))=rVals((minWindow/2));
    gVals(end-(minWindow/2):end)=gVals(end-(minWindow/2));
    rVals(end-(minWindow/2):end)=rVals(end-(minWindow/2));
       fitOptions.StartPoint=[1,-.1,min(rVals)];


f=fit(xVals,rVals,Fexponent,fitOptions);
   fitOptions.StartPoint=[1,-.1,min(gVals)];

g=fit(xVals,gVals,Fexponent,fitOptions);
if f(1)>(max(RvalsAll(i,:))+100)
    f=fit(xVals,rVals,'poly1');
if f.p1>0
    f.p1=0;
end
    
end
if g(1)>(max(GvalsAll(i,:))+100)
    g=fit(xVals,gVals,'poly1');
if g.p1>0
    g.p1=0;
end

  %  g.a=0;
end
%[f.b g.b]
% plot(GvalsAll(i,:))
% hold on
% plot(g)
% imagesc([],[0 1000],ethoTrack')
% alpha(.2);
% ylim([0 max(GvalsAll(i,:))]);
% hold off
% subplot(2,1,2);
% plot(RvalsAll(i,:))
% hold on
% plot(f)
% imagesc([],[0 1000],ethoTrack')
% alpha(.2);
% ylim([0 max(RvalsAll(i,:))]);
% hold off
% drawnow

photoBleachingR(i,:)=f((1:size(RvalsAll,2)))-f(size(RvalsAll,2));

photoBleachingG(i,:)=g((1:size(RvalsAll,2)))-g(size(RvalsAll,2));
    catch
    end
    

end
% photoBleachingR=bsxfun(@minus,photoBleachingR,photoBleachingR(:,end));
% photoBleachingG=bsxfun(@minus,photoBleachingG,photoBleachingG(:,end));
Rvalstemp=RvalsAll-photoBleachingR ;
Rvalstemp(bsxfun(@le,Rvalstemp,quantile(Rvalstemp,.1,2)))=nan;
Rvalstemp=colNanFill(Rvalstemp')';

Gvalstemp=GvalsAll-photoBleachingG ;
Gvalstemp(bsxfun(@le,Gvalstemp,quantile(Gvalstemp,.1,2)))=nan;
Gvalstemp=colNanFill(Gvalstemp')';

runTime=datestr(now,'yyyymmddTHHMMSS');

%%
%f = lsqcurvefit(Fexponent,[.001,-.001,min(photoBleaching)/2],hasPointsTime,photoBleaching')
A=Rvalstemp';

A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);
A(~isnan(A))=A2(~isnan(A));
A(A<-1)=-nan;
R2=A';


A=(Gvalstemp)';
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);
A(~isnan(A))=A2(~isnan(A));
A(A<-1)=-nan;
G2=A';
%A=medfilt2(A,[3,1]);
%A=R2smooth';
% A(R2smooth<.1)=nan;
% A=[A score(1:size(A,1),1:3)];
% A=G2smooth';

A=(Gvalstemp)'./Rvalstemp';
A0=quantile(A,.2,1);
A=bsxfun(@minus, A,A0);
A=bsxfun(@rdivide,A,A0);
%A(~isnan(A))=A2(~isnan(A));
%A(A<-1)=-nan;
A2=colNanFill(A);
A2=medfilt2(A2,[3,1]);
A2=smooth2a(A2,2,0);

Ratio2=A2';

%A=smooth2a(A,0,1);
acorr=corr(A);
atemp=nancov(A)./sqrt(nanvar(A)'*nanvar(A));
acorr(isnan(acorr))=atemp(isnan(acorr));

cg = clustergram(acorr);
cgIdx=str2double(get(cg,'RowLabels'));
[~,cgIdxRev]=sort(cgIdx);
%cgIdx=cgIdx(cgIdx<(max(cgIdx)-2));

%%
fiducialsAll=fiducialPoints(hasPoints);

for i=1:length(fiducialsAll);
        points=fiducialsAll{i};
    points=points(:,1:3);
    emptyX=cellfun(@(x) ~isempty(x),points(:,1));
plotIdx=find(emptyX);
if i==1
    X0=nan(max(plotIdx),length(hasPoints));
Y0=X0;
end
X0(plotIdx,i)=cell2mat(points(plotIdx,1));
Y0(plotIdx,i)=cell2mat(points(plotIdx,2));

end

X0=X0(~rejects,:);
Y0=Y0(~rejects,:);
%%
try
DmatAll=pdistCell(fiducialPoints(hasPoints),find(~rejects),2);
DmatAll=nanmean(DmatAll,3);


corrEdge=-1:.1:1;
dEdge=0:10:370;
[corrD edges mid loc]=histcn([DmatAll(:),acorr(:)],dEdge,-1:.1:1);


corrPlot=accumarray(loc(:,1),acorr(:),size(mid{1}'),@mean);
corrPlotSTD=accumarray(loc(:,1),acorr(:),size(mid{1}'),@std);
corrPlotN=accumarray(loc(:,1),ones(size(acorr(:))),size(mid{1}'))/2;

errorbar(mid{1},corrPlot,corrPlotSTD./sqrt(corrPlotN));
catch
    DmatAll=0;
end


%%
figure
% cmap=flipud(cbrewer('div','RdBu',256,'linear'));
% cmap=interp1([linspace(-.5,0,128) linspace(2/128,2,128)],cmap,linspace(-.5,2,256));

cmap=pmkmp(64,'CubicL');

imagescwithnan(hasPointsTime,1:length(cgIdx),Ratio2(cgIdx,:),...
    cmap,[1 1 1 ]);
caxis([-.3668    2.5])
caxis([-.3668    1.9334])

colorbar off
xlim([0 min(max(hiResFrameTime(hiResV(hiResRange)~=0)),max(hasPointsTime(sum(isnan(G2))<50))')]);
labelSpace=10;
tickLabels=repmat({[]},ceil(size(G2,1)),1);
labelStrings=mat2cell(1:labelSpace:size(G2,1),1,ones(size(1:labelSpace:size(G2,1))));
labelStrings=cellfun(@(x) [num2str(x) ''''] ,labelStrings,'uniformoutput',0);

tickLabels(1:labelSpace:size(G2,1))=labelStrings;

set(gca,'Ytick',1:size(G2,1),'YTickLabel',tickLabels')
xlabel('Time (s)')
ylabel('Neuron ID')
set(gca,'fontSize',15);
savefig([dataFolder filesep runTime 'heatMap' ])

figure
imagesc(hiResFrameTime,1,ethogram');
xlim([0 min(max(hiResFrameTime(hiResV(hiResRange)~=0)),max(hasPointsTime(sum(isnan(G2))<50))')]);

colormap(ethocolormap)
caxis([-1,2])

savefig([dataFolder filesep runTime 'ethogram' ]);



%%
i=40;
figure
plot(R2(cgIdx(i),:),'r')
hold on
plot(G2(cgIdx(i),:),'b')
plot(Ratio2(cgIdx(i),:),'black')
figure
plot(RvalsAll(cgIdx(i),:),'b')
hold on
plot(photoBleachingR(cgIdx(i),:))
hold off


%%

% figure
% hh(1)=subplot(2,1,1);
% imagescwithnan(hasPointsTime,1:length(cgIdx),R2(cgIdx,:),...
%     cmap,[1 1 1 ]);
% hh(2)=subplot(2,1,2);
% imagescwithnan(hasPointsTime,1:length(cgIdx),G2(cgIdx,:),...
%     cmap,[1 1 1 ]);
figure;
imagesc(acorr(cgIdx,cgIdx))

caxis([-1 1]);
axis square
set(gca,'Ytick',1:size(G2,1),'YTickLabel',tickLabels')

colormap(double(flipud(cbrewer('div','RdBu',256,'linear'))))
savefig([dataFolder filesep runTime 'corrMap' ])




%%
save([dataFolder filesep 'heatData'],'G2','R2','Ratio2','ethoTrack','hasPointsTime','acorr','DmatAll','cgIdx','cgIdxRev');

%%
mkdir([dataFolder filesep 'excelData']);

xlswrite([dataFolder filesep 'excelData' filesep 'activity.xlsx'],Ratio2(cgIdx,:)','Activity')
xlswrite([dataFolder filesep 'excelData' filesep 'behavior.xlsx'],ethoTrack,'Behavior')
xlswrite([dataFolder filesep 'excelData' filesep 'time.xlsx'],hasPointsTime,'Time')
xlswrite([dataFolder filesep 'excelData' filesep 'position.xlsx'],[xPosTrack yPosTrack],'Position')

    
%%

progressbar(0)
gflag=1;
bootN=2000;
deltabootGB=zeros(size(G2,1),bootN);
    deltabootGF=deltabootGB;
    deltabootGT=deltabootGB;
    deltabootGP=deltabootGB;
GbootCorrF=deltabootGB;
GbootCorrB=deltabootGB;
GbootCorrT=deltabootGB;
GbootCorrP=deltabootGB;
RbootCorrF=deltabootGB;
RbootCorrB=deltabootGB;
RbootCorrT=deltabootGB;
RbootCorrP=deltabootGB;

    behaviorF=(ethoTrack'==1);
    behaviorP=(ethoTrack'==0);
    behaviorB=(ethoTrack'==-1);   
        behaviorT=(ethoTrack'==2);   


for i=1:size(G2,1);
    progressbar(i/size(G2,1));
    gtemp=smooth(Ratio2(i,:),3)';
    behaviorTemp=ethoTrack';
    
    anans=isnan(gtemp)';
    gtemp=(gtemp(~anans));
    
    behaviorTemp=behaviorTemp(~anans);
    
    behaviorF=(behaviorTemp==1);
    behaviorP=(behaviorTemp==-1|behaviorTemp==2);
    behaviorB=(behaviorTemp==-1);   
        behaviorT=(behaviorTemp==2);   

        bnans=isnan(ethoTrack);
        bnans=bnans(:);
    zTemp=(behaviorZtrack(~(bnans|anans)));
    gtemp2=gtemp(~(bnans|anans));
  %  rtemp2=rtemp(~(bnans|anans));
    
   gzcorr(i)=corr2(gtemp2,zTemp');

  if 0
      
    thetaTemp=behaviorThetaTrack(~(bnans|anans));
    gtemp2=gtemp(~(bnans|anans));

      ggrad=gradient(gtemp2);
   
    [phi(i), thetaCorr(i)]=fminsearch(@(x) -(corr2(ggrad',cos(thetaTemp+x))),0);
   sineTemp=cos(thetaTemp+phi(i));
   
    thetaProj(i)=dot(sineTemp,ggrad)/(sineTemp'*sineTemp);
    
    Ratio3(i,:)=Ratio2(i,:)-sin(behaviorThetaTrack+phi(i))'*thetaProj(i);
  end    
  
    GcorrB(i)=corr2(gtemp,behaviorB);
    GcorrF(i)=corr2(gtemp,behaviorF);
    GcorrP(i)=corr2(gtemp,behaviorP);
    GcorrT(i)=corr2(gtemp,behaviorT);
    deltaGB(i)=nanmean(gtemp(behaviorB))-nanmean(gtemp(~behaviorB));
    deltaGF(i)=nanmean(gtemp(behaviorF))-nanmean(gtemp(~behaviorF));
    deltaGT(i)=nanmean(gtemp(behaviorT))-nanmean(gtemp(~behaviorT));
 %   deltaRP(i)=nanmean(rtemp(behaviorP))-nanmean(rtemp(~behaviorP));
    deltaGP(i)=nanmean(gtemp(behaviorP))-nanmean(gtemp(~behaviorP));
    
    if gflag
            rtemp=smooth(G2(i,:),3)';

    RcorrB(i)=corr2(rtemp,behaviorB);
    RcorrF(i)=corr2(rtemp,behaviorF);
    RcorrP(i)=corr2(rtemp,behaviorP);
    RcorrT(i)=corr2(rtemp,behaviorT);
    deltaRB(i)=nanmean(rtemp(behaviorB))-nanmean(rtemp(~behaviorB));
    deltaRF(i)=nanmean(rtemp(behaviorF))-nanmean(rtemp(~behaviorF));
    deltaRT(i)=nanmean(rtemp(behaviorT))-nanmean(rtemp(~behaviorT));
 %   deltaRP(i)=nanmean(rtemp(behaviorP))-nanmean(rtemp(~behaviorP));
    deltaRP(i)=nanmean(rtemp(behaviorP))-nanmean(rtemp(~behaviorP));
    
    end
    
    
        gTau(i)=min(10,find(diff(autocorr(gtemp-smooth(gtemp,100)',150))>0,1,'first'));

 %   gTau(i)=rTau(i);
    
    %shuffle
        for iTrials=1:bootN
            randOrder=randsamplechunks(length(gtemp),length(gtemp),gTau(i));
        gperm=gtemp(randOrder);
        
    GbootCorrF(i,iTrials)=corr2(gperm,behaviorF(randOrder));
    GbootCorrB(i,iTrials)=corr2(gperm,behaviorB(randOrder));
    GbootCorrP(i,iTrials)=corr2(gperm,behaviorP(randOrder));
    GbootCorrT(i,iTrials)=corr2(gperm,behaviorT(randOrder));
    
    deltabootGB(i,iTrials)=nanmean(gperm(behaviorB))-nanmean(gperm(~behaviorB));
    deltabootGF(i,iTrials)=nanmean(gperm(behaviorF))-nanmean(gperm(~behaviorF));
    deltabootGT(i,iTrials)=nanmean(gperm(behaviorT))-nanmean(gperm(~behaviorT));
    deltabootGP(i,iTrials)=nanmean(gperm(behaviorP))-nanmean(gperm(~behaviorP));
    
      if gflag
    rperm=rtemp((randOrder));

    RbootCorrF(i,iTrials)=corr2(rperm,behaviorF(randOrder));
    RbootCorrB(i,iTrials)=corr2(rperm,behaviorB(randOrder));
    RbootCorrP(i,iTrials)=corr2(rperm,behaviorP(randOrder));
    RbootCorrT(i,iTrials)=corr2(rperm,behaviorT(randOrder));
    
    deltabootRB(i,iTrials)=nanmean(rperm(behaviorB))-nanmean(rperm(~behaviorB));
    deltabootRF(i,iTrials)=nanmean(rperm(behaviorF))-nanmean(rperm(~behaviorF));
    deltabootRT(i,iTrials)=nanmean(rperm(behaviorT))-nanmean(rperm(~behaviorT));
    deltabootRP(i,iTrials)=nanmean(rperm(behaviorP))-nanmean(rperm(~behaviorP));
      end
    
        end
end


    GbootCorrSTDB=quantile(GbootCorrB,.05,2);
    GbootCorrSTDF=quantile(GbootCorrF,.05,2);
    GbootCorrSTDP=quantile(GbootCorrP,.05,2);
    GbootCorrSTDT=quantile(GbootCorrT,.05,2);
          if gflag
    RbootCorrSTDB=quantile(RbootCorrB,.05,2);
    RbootCorrSTDF=quantile(RbootCorrF,.05,2);
    RbootCorrSTDP=quantile(RbootCorrP,.05,2);
    RbootCorrSTDT=quantile(RbootCorrT,.05,2);
          end
          
%     GbootDeltaSTDB(i)=quantile(deltabootGB,.99);
GstdAll=[nanstd(GbootCorrB,[],2) ;nanstd(GbootCorrF,[],2); nanstd(GbootCorrT,[],2)];
RstdAll=[nanstd(RbootCorrB,[],2) ;nanstd(RbootCorrF,[],2); nanstd(RbootCorrT,[],2)];


%%
GRange=range(Ratio2,2);
possibleZ=find(abs(gzcorr)>.5);
deltaCI=.99;
deltaThresh=.2;
nNan=mean(isnan(Ratio2),2)';
nanThresh=.2;
corrThresh=.25;
deltaGBthresh=max(quantile(deltabootGB(:),deltaCI),deltaThresh);

possibleB=find(GbootCorrSTDB'>0 & GcorrB>=corrThresh...
    & deltaGB>deltaGBthresh & nNan<nanThresh);

 [~,ia]=sort(deltaGB(possibleB),'descend');
 possibleB=possibleB(ia);
 possibleB=possibleB(1:min(length(possibleB),16));
% possibleF=find((abs(GcorrF)>1.9*GbootCorrSTDF) & GcorrF>.3 & RGcorr<RGThresh ...
%     & deltaGF>std(0*GbootDeltaSTDF));
deltaGFthresh=max(quantile(deltabootGF(:),deltaCI),deltaThresh);
possibleF=find(GbootCorrSTDF'>0 & GcorrF>corrThresh...
    & deltaGF>deltaGFthresh & nNan<nanThresh);
 [~,ia]=sort(GcorrF(possibleF),'descend');
 possibleF=possibleF(ia);
  possibleF=possibleF(1:min(length(possibleF),6));

  deltaGPthresh=max(quantile(deltabootGP(:),deltaCI),deltaThresh);
possibleP=find(GbootCorrSTDP'>0 & GcorrP>corrThresh  ...
    & deltaGP>deltaGPthresh & nNan<nanThresh);
 [~,ia]=sort(deltaGP(possibleP),'ascend');
 possibleP=possibleP(ia);
 possibleP=possibleP(1:min(length(possibleP),6));
 %possibleP=possibleP(~ismember(possibleP,possibleB));
   deltaGTthresh=max(quantile(deltabootGT(:),deltaCI),deltaThresh);
possibleT=find(GbootCorrSTDT'>0 & GcorrT>corrThresh ...
    & deltaGT>deltaGTthresh & nNan<nanThresh);
 [~,ia]=sort(deltaGT(possibleT),'ascend');
 possibleT=possibleT(ia);
 
 possibleT=possibleT(~ismember(possibleT,possibleB));
 possibleT=possibleT(1:min(length(possibleT),16));

 
 %%
 save([dataFolder filesep runTime 'corrandPossibleCoor2'],'GcorrB','GcorrP','GcorrF','GcorrT',...
 'possibleB','possibleP','possibleF','possibleT','GbootCorrSTDB')
mean([GcorrB(possibleB) GcorrP(possibleT) GcorrF(possibleF)])

%% plot!
spaceVec=[ 1 1 1];
i=-2;
plotOffset=0;
alphaVal=.2;
timeLength=80;%max(hasPointsTime);
ethoSelect=([1;diff(ethogram)]~=0 );
ethoSelect=ethoSelect(~isnan(ethogram));
FplotSelect=~isnan(ethoTrack);
ethoSelect=smooth(ethoSelect,15)>0;
ethoSelect(end)=true;
    ethoHeight=10;
%%
topAll=10;
for class=1:3
    space=spaceVec(class);
    plotOffset=i+2+plotOffset;
    plotOffset=-.5;
    switch class
        case 1
select=[possibleF ];
        case 2
            select=possibleT;
        case 3
            select=possibleB;
    end
    space=(space*topAll/(length(select)+2));
    if ~isempty(select);
%select=select(1:min(6,length(select)));

for i=1:length(select)
        subplot(3,1,class);
        plot(frameTimeTrack(FplotSelect),medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset),'black','linew',1)
hold on
topPoint=max(medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset));
topAll=max(topAll,topPoint);
%plot(R2(possibleB(i),:)'+space*i,'r')

end

    subplot(3,1,class);
%ylim([0 ethoHeight])
set(gca,'YTick',1,'fontsize',12)
xlabel('Time(s)')
ylabel('\Delta F /F0')
  text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
        ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
        'fontsize',13);
    xlim([0 max(frameTimeTrack(FplotSelect)+15)])
    %axis square
    
    h=imagesc(hiResFrameTime,[1 10],ethogram');
colormap(ethocolormap)
 h.Face.ColorType='truecoloralpha';
set(h,'alphaData',ones(size(ethogram'))*alphaVal)
ylim([0,topAll]);
% 

% h(1)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)<0),'facecolor',bcolor);
% h(2)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==0),'facecolor',pausecolor);
% h(3)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==1),'facecolor',fcolor);
% h(4)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==2),'facecolor',turncolor);
% drawnow
% for ii=1:length(h)
%     
%     subplot(3,1,class);
%     h(ii).EdgeColor='none';
%     h(ii).Face.ColorType='truecoloralpha';
%     h(ii).Face.ColorData(4)=alphaVal*255;
% end
% 
% set(gca,'YTick',space,'fontsize',12)
% xlabel('Time(s)')
% ylabel('\Delta F /F0')
%   text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
%         ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
%         'fontsize',13);
%     xlim([0 timeLength])
% h(1)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)<0),'facecolor',bcolor);
% h(2)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==0),'facecolor',pausecolor);
% h(3)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==1),'facecolor',fcolor);
% h(4)=area(hiResFrameTime(ethoSelect),ethoHeight*(ethogram(ethoSelect)==4),'facecolor',turncolor);
% drawnow
%     axis square
% for i=1:length(h)
%     
% 
%     h(i).EdgeColor='none';
%     h(i).Face.ColorType='truecoloralpha';
%     h(i).Face.ColorData(4)=alphaVal*255;
% end

    end
end


%savefig([dataFolder filesep runTime 'plots' ])


  
%%
% masterData=load([dataFolder filesep 'refStackStraight']);
% fiducials=masterData.masterFiducials;

nPoints=cellfun(@(x) length(cell2mat(x)) ,fiducialPoints);
allPoints=find(nPoints==max(nPoints));
fiducials=(fiducialPoints{allPoints(1)});
fiducials=fiducials(~rejects,[1 2 3]);
fiducials=cell2mat(fiducials);

figure
colorBalls=.3+zeros(size(fiducials));
colorBalls((possibleT),3)=1;
colorBalls((possibleB),1)=1;
colorBalls((possibleF),2)=1;
transp=ones(1,length(fiducials));
transp=transp*.3;
transp(([possibleT possibleB possibleF]))=1;
scatter3sph(fiducials(:,1),fiducials(:,2),60*fiducials(:,3),...
    'size',7,'color',colorBalls,'transp',transp);axis equal off

savefig([dataFolder filesep runTime 'balls' ])

%% check out track of moving fiducial postions


hiResS2low=load('Y:\CommunalCode\3dbrain\registration\20141212HighResS2LowResFluorBeads.mat');
lowRes=load('Y:\CommunalCode\3dbrain\registration\20141212LowResBehavior2Fluor.mat');

%%
for i=1:2%:length(hasPoints);

    fiducialsI=fiducialPoints{hasPoints(i)};
    plotIdx=find(cell2mat(cellfun(@(x) ~isempty(x),fiducialsI(:,1),'uniform',0)));
    fiducialsI=cell2mat(fiducialsI);
  fiducialsI(:,1)=fiducialsI(:,1);
  fiducialsI(:,2)=fiducialsI(:,2);
  if i==1
M=mean(fiducialsI);
  end
X=-(fiducialsI(:,2)-M(2))/4000+yPosStage(fiducialsI(:,end));
    Y=(fiducialsI(:,1)-M(1))/4000-xPosStage(fiducialsI(:,end));
    
    
    bfFrame=round(mean(bfIdxLookup(fiducialsI(:,end))));
    CLFrame=bfFrame-offset;
    CL=centerline(:,:,CLFrame);
    [CL2(:,2),CL2(:,1)]=transformPointsInverse(lowRes.t_concord,CL(:,2),CL(:,1));

[CL2(:,1),CL2(:,2)]=transformPointsForward(hiResS2low.t_concord,CL2(:,2),CL2(:,1));
% CL2(:,1)=CL2(:,1)-1200;
 CL2(:,2)=CL2(:,2)-600;

scatter3(fiducialsI(:,1),fiducialsI(:,2),30*fiducialsI(:,3))
hold on
%plot(CL2(:,1),CL2(:,2));


    CLX=-(CL2(:,2)-M(2))/4000+yPosStage(fiducialsI(10,end));
    CLY=(CL2(:,1)-M(1))/4000-xPosStage(fiducialsI(10,end));
%     if i==1
% plot(yPosStage,-xPosStage);
% hold on
%     h=scatter(X,Y,'x');
% %text(X,Y,cellstr(num2str(plotIdx)))
% h2=plot(CLX,CLY);
%     else
%          h.XData=X;
%         h.YData=Y;
% %         
%         h2.XData=CLX;
%         h2.YData=CLY;
%     end
        

axis equal
%     xlim([-5 ,0]);
%      ylim([-5,1]);
    pause(.3)
   

end


%%
for i=1:length(hasPoints);

    fiducialsI=fiducialPoints{hasPoints(i)};
    plotIdx=find(cell2mat(cellfun(@(x) ~isempty(x),fiducialsI(:,1),'uniform',0)));
    fiducialsI=cell2mat(fiducialsI);
  fiducialsI(:,1)=fiducialsI(:,1);
  fiducialsI(:,2)=fiducialsI(:,2);
  if i==1
M=mean(fiducialsI);
  end
X=-(fiducialsI(:,2)-M(2))/4000+yPosStage(fiducialsI(:,end));
    Y=(fiducialsI(:,1)-M(1))/4000-xPosStage(fiducialsI(:,end));
    
    
    bfFrame=round(mean(bfIdxLookup(fiducialsI(:,end))));
    CLFrame=bfFrame-offset;
    CL=centerline(:,:,CLFrame);
    [CL2(:,2),CL2(:,1)]=transformPointsInverse(lowRes.t_concord,CL(:,2),CL(:,1));

[CL2(:,1),CL2(:,2)]=transformPointsForward(hiResS2low.t_concord,CL2(:,2),CL2(:,1));
% CL2(:,1)=CL2(:,1)-1200;
 CL2(:,2)=CL2(:,2)-600;
 
 
%  
%  
%  
% 
% scatter3(fiducialsI(:,1),fiducialsI(:,2),30*fiducialsI(:,3))
% hold on
%plot(CL2(:,1),CL2(:,2));


    CLX=-(CL2(:,2)-M(2))/4000+yPosStage(fiducialsI(10,end));
    CLY=(CL2(:,1)-M(1))/4000-xPosStage(fiducialsI(10,end));
    
    
   [t_vector,b_vector,n_vector]=tbnVector([CLX,CLY,0*CLX]); 
  
   Dmat=pdist2([CLX,CLY],[X,Y]);
   [~,Dmin]=min(Dmat);
   S=dot(n_vector(Dmin,1:2),([X,Y]-[CLX(Dmin),CLY(Dmin)]),2);
   scatter(Dmin,S)
   
%     if i==1
% plot(yPosStage,-xPosStage);
% hold on
%     h=scatter(X,Y,'x');
% %text(X,Y,cellstr(num2str(plotIdx)))
% h2=plot(CLX,CLY);
%     else
%          h.XData=X;
%         h.YData=Y;
% %         
%         h2.XData=CLX;
%         h2.YData=CLY;
%     end
        

%     xlim([-5 ,0]);
      ylim([-.03,.03]);
    pause(.3)
   

end




