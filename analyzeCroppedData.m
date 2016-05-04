%dataFolder=uipickfiles;
%dataFolder=dataFolder{1};
dataFolder=uipickfiles();
dataFolder=dataFolder{1};
%'F:\20141212\BrainScanner20141212_145951\';
runTime=datestr(now,'yyyymmddTHHMMSS');


%% load Fiducials file
 [fiducialPoints,z2ImageIdxOffset]=loadFiducialPoints(dataFolder);
 %%
 try
    alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
catch
display('Select Low Res Alignment')

lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
%lowResFluor2BF=load('Y:\CommunalCode\3dbrain\registration\20141212LowResBehavior2Fluor.mat');
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);


display('Select Hi to Low Fluor Res Alignment')
Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load(Hi2LowResF{1});
%Hi2LowResF=load('Y:\CommunalCode\3dbrain\registration\20141212HighResS2LowResFluorBeads.mat');


% display('Select Hi to Low Res Alignment')
% 
% Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
% Hi2LowRes=load(Hi2LowRes{1});
% t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');
 display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
%S2AHiRes=load('Y:\CommunalCode\3dbrain\registration\20141212HiResS2A.mat');
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;

%%


%%
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;

save([dataFolder filesep 'alignments'],'alignments');

end

%%
imSize=[1200,600];
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
[bfTime,hiResFrameTime,hasPoints,bfRange,hiResRange,hasPointsTime,lookup]=dataTimeAlignment(dataFolder,fiducialPoints);
 
nTimes=length(hasPoints);

%% load centerline data
[centerline, offset,eigenProj, CLV,wormCentered]=loadCLBehavior(dataFolder);


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



hiResCLposition=[interp1(1:length(centerLinePosition),centerLinePosition(:,1),lookup.BF2hiRes,'nearest','extrap'),...
    interp1(1:length(centerLinePosition),centerLinePosition(:,2),lookup.BF2hiRes,'nearest','extrap')];
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


behaviorZtrack=behaviorZ(bfRange);
projectionsTrack=projections(bfRange,:);
behaviorThetaTrack=behaviorTheta(bfRange);



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

filterKernal2=gausswin(500);
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

hiResV=v;%gradient(D,50)./median(gradient(hiResData.frameTime));
hiResV=hiResV.*hiResCLbehavior;
hiResV=smooth(hiResV,500);
%%
vVectorSmooth=[gradient(smooth(xPos,1000)) gradient(smooth(yPos,1000))];
vAngle=(atan2(vVectorSmooth(:,2),vVectorSmooth(:,1)));
vAngleTrack=vAngle((find(diff(hiResData.stackIdx)>0)));
vAngleTrack=vAngleTrack(hasPoints);
[vHist,vEdges,vLocs] = histcounts(vAngleTrack,-pi:pi/12:pi);
polar(vEdges(2:end),vHist)

if 1
for iNeuron=1:size(Ratio2,1);
angleCorr(iNeuron,:)=accumarray(vLocs,Ratio2(iNeuron,hasPoints),[],@nanmean);

end
angleMean=mean(angleCorr);
angleSTD=std(angleCorr);
angleResidual=bsxfun(@minus,angleCorr,angleMean);
angleResidual=bsxfun(@rdivide,angleResidual,angleMean);

angleResidual=(sum(angleResidual.^2,2));
end


% 
% figure;
% plot(hiResData.frameTime(hiResRange)-hiResData.frameTime(hiResRange(1)),hiResV(hiResRange),...
%     'color',[106 61 154]/256,'linew',2);
% axis tight;
% ylim([-.15, .15])
% xlabel('Time (s)');
% ylabel('Velocity (mm/s)');
% set(gca,'box','off');
% %xaxescenter
% set(gca,'fontSize',15);


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
ethoZ=hiResBehaviorZ(hiResRange);
ethoZ=ethoZ-nanmean(ethoZ);
ethoZ(isnan(ethoZ))=0;
Vcluster=smooth(hiResV(hiResRange),100);
idx=sign(Vcluster);
idx(abs(Vcluster)<.00005)=0;
ethogram=idx;
pausing=find(ethogram==0);

idxgmm=kmeans(ethoZ,2,'start',[ 0 5]');
%%
ethogram((((abs(ethoZ)>2*std(ethoZ))& (ethogram>=0)))|abs(ethoZ)>10)=2;
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
centerlineShowFrame=bfRange(round(linspace(1,length(bfRange),5)));
centerlineShowFrame(centerlineShowFrame>size(centerline,3))=size(centerline,3);
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
    axis equal tight off

%%
%imagesc(hiResFrameTime,1,ethogram');

% h=colorplot(hiResFrameTime/max(hiResFrameTime)*2500,1,ethogram);
% colormap(ethocolormap)

%h=colorplot(hiResFrameTime/max(hiResFrameTime)*2500,1,hiResFrameTime);

%set(h,'linew',50);

    

heatMapGeneration_spaceCorr(dataFolder,fiducialPoints,alignments,runTime)
load([dataFolder filesep 'heatData' runTime]);

nanMap=isnan(rRaw)|isnan(gRaw)| rRaw<40;
nanMap2=imopen(imclose(nanMap,[1 1 1 1 1]),[1 1 1 1 1]);
Ratio2(nanMap2)=nan;
R2(nanMap2)=nan;
G2(nanMap2)=nan;
%%
if 0
    [FX,FY,FZ]=fiducials2mat(fiducialPoints,size(Ratio2,1));
%%
    hg=histcn([FX(:) FY(:)],'AccumData',G2(:),'FUN',@(x) trimmean(x,20));
    hr=histcn([FX(:) FY(:)],'AccumData',R2(:),'FUN',@(x) trimmean(x,20));
    hratio=histcn([FX(:) FY(:)],'AccumData',Ratio2(:),'FUN',@(x) trimmean(x,20));
hg(hg==0)=nan;
hr(hr==0)=nan;
hratio(hratio==0)=nan;
    
h(1)=subplot(1,3,1);
[h hcb]=imagescwithnan(hr,parula,[1 1 1]);
delete(hcb)
axis off square
caxis([-.5 1])
h(2)=subplot(1,3,2);
[h hcb]=imagescwithnan(hg,parula,[1 1 1]);
delete(hcb)
axis off square
caxis([-.5 1])
h(3)=subplot(1,3,3);
[h hcb]=imagescwithnan(hratio,parula,[1 1 1]);
delete(hcb)
axis off square
caxis([-.5 1])

%hlink = linkprop(h,{'CameraPosition','CameraUpVector','xlim','ylim'});
 %   view(0,90)


good=~isnan(Ratio2(:)+FY(:)) & FX(:)>10 & G2(:)<10;
fs=fit([FX(good),FY(good)],Ratio2(good),'poly33');
RatioCorrection=fs(FX,FY);
RatioCorrection=bsxfun(@rdivide,RatioCorrection,nanmean(RatioCorrection,2));
%%
Ratio3=Ratio2./RatioCorrection;

%%
load([dataFolder filesep 'pointStatsNew']);
for iRef=100:length(pointStatsNew)
RefPS=pointStatsNew(iRef);
RefPSPresent=~isnan(RefPS.trackIdx);
trackIdx=RefPS.trackIdx(RefPSPresent);
straightPoints=RefPS.straightPoints(RefPSPresent,:);
[~,ia]=sort(trackIdx);
straightPoints=straightPoints(ia,:);
scatter3(straightPoints(:,1),straightPoints(:,2),straightPoints(:,3),[],cgIdx,'fill');axis equal off tight
drawnow
end


%%
iRef=200;
RefPS=pointStatsNew(iRef);
RefPSPresent=~isnan(RefPS.trackIdx);
trackIdx=RefPS.trackIdx(RefPSPresent);
straightPoints=RefPS.straightPoints(RefPSPresent,:);
[~,ia]=sort(trackIdx);
straightPoints=straightPoints(ia,:);

else
Ratio3=Ratio2;
end

save([dataFolder filesep 'heatData'],'G2','R2','gRaw','rRaw','Ratio2',...
    'gPhotoCorr','rPhotoCorr','acorr','cgIdx','cgIdxRev','ethoTrack','hasPointsTime');

%%
figure
% cmap=flipud(cbrewer('div','RdBu',256,'linear'));
% cmap=interp1([linspace(-.5,0,128) linspace(2/128,2,128)],cmap,linspace(-.5,2,256));
cRange=[-.25 1.5];
cmap=pmkmp(64,'CubicL');
Ratio3plot=Ratio3;
Ratio3plot(Ratio3<cRange(1))=cRange(1);
Ratio3plot(Ratio3> cRange(2))= cRange(2);
h(1)=imagescwithnan(hasPointsTime,1:length(cgIdx),Ratio3plot(cgIdx,:),...
    cmap,[1 1 1 ]);
% caxis([-.3668    2.5])
% caxis([-.3668    1.9334])

colorbar off
xlim([0 max(hasPointsTime)]);
labelSpace=10;
tickLabels=repmat({[]},ceil(size(G2,1)),1);
labelStrings=mat2cell(1:labelSpace:size(G2,1),1,ones(size(1:labelSpace:size(G2,1))));
labelStrings=cellfun(@(x) [num2str(x) ''''] ,labelStrings,'uniformoutput',0);

tickLabels(1:labelSpace:size(G2,1))=labelStrings;

set(gca,'Ytick',1:size(G2,1),'YTickLabel',tickLabels')
xlabel('Time (s)')
ylabel('Neuron ID')
set(gca,'fontSize',15);
savefig([dataFolder filesep 'heatMap' ])

figure
h(2)=imagesc(hiResFrameTime,1,ethogram');
set(gca,'fontSize',15);
set(gca,'YTick','')
%xlim([0 min(max(hiResFrameTime(hiResV(hiResRange)~=0)),max(hasPointsTime(sum(isnan(G2))<50))')]);

colormap(ethocolormap)
caxis([-1,2])
xlim([0 max(hasPointsTime)]);
linkaxes([h(1).Parent h(2).Parent],'x')

savefig([dataFolder filesep 'ethogram' ]);



%%

figure
hh(1)=subplot(2,1,1);
imagescwithnan(hasPointsTime,1:length(cgIdx),R2(cgIdx,:),...
    cmap,[1 1 1 ]);
hh(2)=subplot(2,1,2);
imagescwithnan(hasPointsTime,1:length(cgIdx),G2(cgIdx,:),...
    cmap,[1 1 1 ]);
%%
figure;
imagesc(acorr(cgIdx,cgIdx))

caxis([-1 1]);
axis square
set(gca,'Ytick',1:size(G2,1),'YTickLabel',tickLabels')

colormap(double(flipud(cbrewer('div','RdBu',256,'linear'))))
savefig([dataFolder filesep runTime 'corrMap' ])




%%
save([dataFolder filesep 'heatData'],'G2','R2','Ratio2','ethoTrack','hasPointsTime','acorr','cgIdx','cgIdxRev','gRaw','rRaw');

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
    
   gzcorr(i)=corr2(gtemp2,zTemp);

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
  possibleF=possibleF(1:min(length(possibleF),12));

  deltaGPthresh=max(quantile(deltabootGP(:),deltaCI),deltaThresh);
possibleP=find(GbootCorrSTDP'>0 & GcorrP>corrThresh  ...
    & deltaGP>deltaGPthresh & nNan<nanThresh);
 [~,ia]=sort(deltaGP(possibleP),'descend');
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
close all
subplotFlag=0;
topAll=30;
nPlots=12;
for class=1:3
    space=spaceVec(class);
    plotOffset=i+2+plotOffset;
    plotOffset=-.5;
    switch class
        case 1
select=[possibleP ];
        case 2
            select=possibleT(1:min(nPlots,length(possibleT)));
        case 3
            select=possibleB(1:min(nPlots,length(possibleB)));
    end
    space=(space*topAll/(length(select)+1))
    if ~isempty(select);
%select=select(1:min(6,length(select)));
clear h2
if subplotFlag
  h2(class)=      subplot(3,1,class);
else
    figure(class)
end
for i=1:length(select)
        plot(frameTimeTrack(FplotSelect),medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset),'black','linew',1)
hold on
topPoint=max(medfilt1(Ratio2(select(i),FplotSelect),5)'+space*(i+plotOffset));
%plot(R2(possibleB(i),:)'+space*i,'r')

end

%    subplot(3,1,class);
%ylim([0 ethoHeight])
set(gca,'YTick',1,'fontsize',12)
xlabel('Time(s)')
ylabel('\Delta R /R0')
  text( max(frameTimeTrack(FplotSelect)+5)*ones(size(select)), space*plotOffset+[space:space:space*(i)],cellstr(num2str(cgIdxRev(select))),'VerticalAlignment'...
        ,'middle', 'HorizontalAlignment','left','color',[0 0 0],...
        'fontsize',13);
    xlim([0 max(frameTimeTrack(FplotSelect)+15)])
    %axis square
    
    h=imagesc(hiResFrameTime,[1 topPoint],ethogram');
colormap(ethocolormap)
 h.Face.ColorType='truecoloralpha';
set(h,'alphaData',ones(size(ethogram'))*alphaVal)
ylim([0,topAll]);
% 


    end
end


%savefig([dataFolder filesep runTime 'plots' ])


  
%%
% masterData=load([dataFolder filesep 'refStackStraight']);
% fiducials=masterData.masterFiducials;

nPoints=cellfun(@(x) length(cell2mat(x)) ,fiducialPoints);
allPoints=find(nPoints==max(nPoints));
fiducials=straightPoints;
figure
colorBalls=.2+zeros(size(fiducials));
b=possibleP%(1:min(nPlots,length(possibleP)))
t=[];possibleT(1:min(nPlots,length(possibleT)));
f=[];possibleF(1:min(nPlots,length(possibleF)));
colorBalls((t),:)=repmat([55,126,184]/256,length(t),1);
colorBalls((b),:)=repmat([228,26,28]/256,length(b),1);
colorBalls((f),:)=1;
transp=ones(1,length(fiducials));
transp=transp*.2;
transp([t b f])=1;
scatter3sph(fiducials(:,1),fiducials(:,2),fiducials(:,3),...
    'size',6,'color',colorBalls,'transp',transp);axis equal off
text(fiducials(b,1),fiducials(b,2)+10*[-3 1 1 1 1]',fiducials(b,3)+10,...
    cellstr(num2str(cgIdxRev(b))),'FontSize',18)

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

