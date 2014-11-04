%%  load syncing data
dataFolder=uipickfiles('filterspec','E:\DATA\3DwormData\');%'E:\DATA\3DwormData\BrainScanner20140911_182843';
dataFolder=dataFolder{1};
hiResData=highResTimeTraceAnalysis(dataFolder,600,1200);


hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
hi2lowTimeOffset=hiResFlashTime-bfFlashTime;
hiResData.frameTime=hiResData.frameTime-hi2lowTimeOffset(1);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

fluor2bfTimeOffset=fluorFlashTime-bfFlashTime;
fluorAll.frameTime=fluorAll.frameTime-fluor2bfTimeOffset(1);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

%works amazingly well, flash frames off by less than .05 seconds. 


%% get Z timing
Zall=hiResData.Z;
timeAll=hiResData.frameTime;
changes=[0;diff(Zall)<-.6];
zSpikes=[0; diff(changes)<-.5];
stackIdx=cumsum(zSpikes);
%% load centerline data
centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile={centerLineFile.name}';
if length(centerLineFile)>1
    centerlineFile=uipickfiles('FilterSpec',dataFolder);
    centerline=load(centerlineFile{1},'centerline');
    centerline=centerline.centerline;
else
centerline=load([dataFolder filesep centerLineFile],'centerline');
centerline=centerline.centerline;
end
%% load alignment data
display('Select Low Res Alignment')
lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);
display('Select Hi to Low Res Alignment')

Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowRes=load(Hi2LowRes{1});
t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');
display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;

%% load background image
try
backgroundImageFile=dir([dataFolder filesep '*backgroundImage*']);
if isempty(backgroundImageFile) || length(backgroundImageFile)>1
    display('Select Background Mat File')
    backgroundImageFile=uipickfiles('FilterSpec',dataFolder);
    backgroundImageFile=backgroundImageFile{1};
else
    backgroundImageFile=[dataFolder filesep backgroundImageFile.name];
end
backgroundImage=load(backgroundImageFile);
backgroundImage=backgroundImage.backgroundImage;
catch
    backgroundImage=0;
end




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

Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

%% make lookup tables for indices
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));

%%
stretchSize=25;
%lowResFolder=[dataFolder filesep 'lowResFolder'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder3D'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder3D'];
metaFolder=[dataFolder filesep  'metaDataFolder3D'];
mkdir(metaFolder);
mkdir(hiResActivityFolder);
mkdir(hiResSegmentFolder);
frames=1750:length(hiResData.frameTime); %start 1750, 12000 good too, 13000 for 3d
frames(ismember(frames,hiResData.flashLoc))=[];
movieFlag=0;
plotFlag=0;
saveFlag=1;
imH=NaN(1,4);
lineH=NaN(1,4);

hiResCorrection=[0,0];
spikeBuffer=5;
meanSliceFiltLevel=4;
yWindow=300;
xWindow=200;
%%
%progressbar(0,0)
stackList=2:max(stackIdx);
overwriteFlag=1;
if ~overwriteFlag;
    d=dir([hiResActivityFolder filesep '*.tif']);
    d={d.name}';
    oldIdx=cell2mat(cellfun(@(x) str2double(x(6:10)),d,'UniformOutput',false));
    stackList=stackList((~ismember(stackList,oldIdx)'));
end


errors=struct('me',[]);
errors=repmat(errors,max(stackIdx)*2,1);
metaData=struct('time',[]);
metaData=repmat(metaData,max(stackIdx)*2,1);
%%
parfor i=1:length(stackList)
    tic
    iStack=stackList(i);
    display([num2str(iStack),';' num2str(i)]);
    %%
    try
        imageIdx=find(stackIdx==iStack);
    imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
worm=[];
zPos=hiResData.Z(imageIdx);
time=hiResData.frameTime(imageIdx);
activity=[];
highResInterpStack2=[];
highResInterpStack=[];
hiResMean=[];
newYHiFout=[];
newXHiFout=[];
segmentStack=[];
activityStack=[];

    %%
    Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

    
    for iSlice=1:length(imageIdx)
      %  progressbar(iStack/max(stackIdx),iSlice/length(imageIdx))
        iFrame=imageIdx(iSlice);
    iTime=hiResData.frameTime(iFrame);
    %interpolate using time to get low res idx 
    bfIdx=bfIdxLookup(iFrame);
fluorIdx=fluorIdxLookup(iFrame);
hiResIdx=hiResLookup(iFrame);

  status=fseek(Fid,2*hiResIdx*1024^2,-1);

  pixelValues=fread(Fid,1024^2,'uint16',0,'l');

hiResImage=(reshape(pixelValues,1024,1024)');
hiResImage=hiResImage-backgroundImage;
hiResImage(hiResImage<0)=0;

 activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
 segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
 activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
activityChannel=pedistalSubtract(activityChannel,5);
segmentStack(:,:,iSlice)=segmentChannel;
activityStack(:,:,iSlice)=activityChannel;
% 
%  hiResImage2=imwarp(hiResImage,Hi2LowRes.t_concord,...
%     'OutputView',Hi2LowRes.Rsegment);

fluorFrame=read(fluorVidObj,round(fluorIdx));
bfFrame = read(behaviorVidObj,round(bfIdx));
fluorFrame=fluorFrame(:,:,1);
bfFrame=bfFrame(:,:,1);
fluorFrame=pedistalSubtract(fluorFrame);

% fluorFrame2=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
%     'OutputView',lowResFluor2BF.Rsegment);
% fluorFrame2=imwarp(fluorFrame2,Hi2LowRes.t_concord,  'OutputView',Hi2LowRes.Rsegment);
% fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

%  bfFrame=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
%     'OutputView',lowResFluor2BF.Rsegment);
%%
CLcurrent=[interp2(squeeze(centerline(:,1,:))',1:100,repmat(bfIdx,1,100))',...
interp2(squeeze(centerline(:,2,:))',1:100,repmat(bfIdx,1,100))'];
% CLcurrent=[csaps(1:100,CLcurrent(:,1),p,1:100,W)'...
%     csaps(1:100,CLcurrent(:,2),p,1:100,W)'];
CLcurrent=[interp1(CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')',...
    interp1(CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')'];

endpts=CLcurrent([1,length(CLcurrent)],:);
[maxPointx,maxPointy]=find(fluorFrame==max(fluorFrame(:)),1,'first');
fluorFrameThresh=im2bw(fluorFrame,graythresh(fluorFrame));
fluorStats=regionprops(fluorFrameThresh,fluorFrame,'Area','WeightedCentroid');
fluorAreas=[fluorStats.Area];
maxPoint=fluorStats(find(fluorAreas==max(fluorAreas),1,'first')).WeightedCentroid;

%maxPoint=[maxPointx,maxPointy];
tip2maxDistance=pdist2(endpts,maxPoint);
% if tip2maxDistance(2)<tip2maxDistance(1)
%     CLcurrent=flipud(CLcurrent);
% end


CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
totLengthPix=max(CLdistance);
CLcurrent=[interp1(CLdistance,CLcurrent(:,1),1:1:totLengthPix)',...
    interp1(CLdistance,CLcurrent(:,2),1:1:totLengthPix)'];



[~,newX,newY]=wormStraightening(CLcurrent,[],80,10);

    [~,newXHi,newYHi]=wormStraightening(CLcurrent(1:min(Inf,length(CLcurrent)),:),[],80,10);
[Xgrid, Ygrid]=meshgrid(1:.02:size(newXHi,2),1:.02:size(newXHi,1));
newXHi=interp2(newXHi,Xgrid,Ygrid);
newYHi=interp2(newYHi,Xgrid,Ygrid);    
    
% CLnew=([fliplr(CLcurrent) ones(length(CLcurrent),1)]*lowResBF2FluorT.T);
% CLnew=CLnew(:,[2,1]);

   [CLnewY,CLnewX]=transformPointsInverse(lowResFluor2BF.t_concord...
       , CLcurrent(:,2),CLcurrent(:,1));
   
   CLnew=[CLnewX CLnewY];
   [CLHighResY,CLHighResX]=transformPointsForward(Hi2LowRes.t_concord...
       , CLcurrent(:,2),CLcurrent(:,1));
   
   CLHighRes=[CLHighResX, CLHighResY];
   CLHighRes=bsxfun(@plus,CLHighRes,hiResCorrection);
   CLHighRes=CLHighRes(CLHighResX<1024 & CLHighResY<1024 & ...
       CLHighResX>0 & CLHighResY>512,:);
   
   [ newYHiF,newXHiF]=transformPointsForward(Hi2LowRes.t_concord...
       , newYHi,newXHi);
   newYHiF=newYHiF-rect1(1)-1;
   newXHiF=newXHiF-rect1(2)-1;
%    n
%       [ newYHiA,newXHiA]=transformPointsInverse(S2AHiRes.t_concord...
%        , newYHi,newXHi);
   %%
   highResInterp=interp2(segmentChannel,newYHiF,newXHiF);
   highResActivityInterp=interp2(activityChannel,newYHiF,newXHiF);

     [ newYHiFlo,newXHiFlo]=transformPointsInverse(lowResFluor2BF.t_concord...
       , newYHi,newXHi);
   
   lowResFluorHiInterp=interp2(double(fluorFrame),newYHiFlo,newXHiFlo);
      lowResFluorHiInterp=pedistalSubtract(lowResFluorHiInterp,5);
      lowResFluorHiInterp(isnan(lowResFluorHiInterp))=0;

    
    
  fluorThresh=lowResFluorHiInterp>(max(lowResFluorHiInterp(:))/2);
  if any(fluorThresh(:))
stats=regionprops(fluorThresh,lowResFluorHiInterp,'Area','WeightedCentroid');
centroid=stats([stats.Area]==max([stats.Area])).WeightedCentroid;
 maxX=round(centroid(1));
 maxY=round(centroid(2));
  else
  maxX=size(lowResFluorHiInterp,2)/2;
  maxY=size(lowResFluorHiInterp,1)/2;
      
  end
% maxY=find(sum(lowResFluorHiInterp,2)==max(sum(lowResFluorHiInterp,2)));
%  maxX=find(sum(lowResFluorHiInterp,1)==max(sum(lowResFluorHiInterp,1)));

lowResFluorHiInterp=rectCrop(lowResFluorHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
highResInterp=rectCrop(highResInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
highResActivityInterp=rectCrop(highResActivityInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
newYHiFout(:,:,iSlice)=rectCrop(newYHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
newXHiFout(:,:,iSlice)=rectCrop(newXHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);

highResInterp=pedistalSubtract(highResInterp);
%lowResFluorHiInterpStack(:,:,iSlice)=pedistalSubtract(lowResFluorHiInterp);
highResInterpStack(:,:,iSlice)=highResInterp;
%highResActivityInterpStack(:,:,iSlice)=pedistalSubtract(highResActivityInterp);
hiResMean(iSlice)=nanmean((highResInterp(:)));

    end
    fclose(Fid);
    %%
[~,valleyPos]=findpeaks(-hiResMean,'THRESHOLD',.5);
goodSlice=true(1,iSlice);
goodSlice(valleyPos)=false;
 highResInterpStack=highResInterpStack(:,:,goodSlice);
 
 activityStack=activityStack(:,:,goodSlice);
 segmentStack=segmentStack(:,:,goodSlice);
 newXHiFout=newXHiFout(:,:,goodSlice);
  newYHiFout=newYHiFout(:,:,goodSlice);
 
 zPos=zPos(goodSlice);
 time=time(goodSlice);
imageIdx=imageIdx(goodSlice);
  highResInterpStack(isnan(highResInterpStack))=0;
  
  
stackSize=size(highResInterpStack,3);

 %%
 subWorm=highResInterpStack(1:size(highResInterpStack,1)/2,:,:);
subWorm2=highResInterpStack(size(highResInterpStack,1)/2:end,:,:);
noseProfile=smooth(squeeze(sum(sum(subWorm(:,:,:),1),2)),3);
  [ymax,imax,ymin,imin]=extrema(noseProfile);
 if length(imax)>1
  maxPeaksPos=sort(imax(1:2));
midPlane=min(imin(imin>maxPeaksPos(1) & imin<maxPeaksPos(2)));
 else
     midPlane=0;
 end
 
if midPlane>stackSize*.6 || midPlane<stackSize*.4
     sliceIntensity=(squeeze(sum(sum(subWorm2(20:end,:,:),1),2)));
[~,midPlane]=findpeaks(sliceIntensity);
tmp=abs(midPlane-stackSize/2);
     midPlane=(midPlane(tmp==min(tmp)));
     midPlane=max(midPlane)+1;
    
end
%%
subBotWorm=subWorm(:,:,1:midPlane);
subTopWorm=subWorm(:,:,midPlane:end);
botWorm=highResInterpStack(:,:,1:midPlane);
subBotWormProj=max(subBotWorm,[],3);
botWormProj=max(botWorm,[],3);
subTopWormProj=max(subTopWorm,[],3);


Pbot=WormBrain3Points(subBotWormProj);
Ptop=WormBrain3Points(subTopWormProj);

if length(Pbot)>length(Ptop)
    Ptop=[Pbot(1,:);Ptop];
elseif length(Ptop)>length(Pbot)
    Pbot=[Ptop(1,:);Pbot];

end



%%
subWorm2Proj=normalizeRange(sum(subWorm2(:,:,midPlane+(-2:2)),3));

  [ymax,imax,ymin,imin]=extrema(smooth(squeeze(sum(subWorm2Proj,2)),5));
brainSplit=imin(ymin<max(ymax)/2);
brainSplit=min(brainSplit);

subWorm2Proj=subWorm2Proj(brainSplit:end,:);
%subWorm2Proj=smooth2a(subWorm2Proj,10,10);

subBW=im2bw(subWorm2Proj,graythresh(subWorm2Proj));
subBW=AreaFilter(subBW,10);
subBW(:,250:end)=0;
[y,x]=find(subBW);
[y,~,ib]=unique(y);
x=accumarray(ib,x,[],@min);
x=x(20:end);
y=y(20:end);
y=y+brainSplit+size(highResInterpStack,1)/2;
y(x>xWindow)=[];
x(x>xWindow)=[];


x=smooth(x,5);
y=smooth(y,5);
%%
  [ymax,imax,ymin,imin]=extrema(x);

 imin=sort(imin);
 imin=imin(2:end);
 
  x=x([imin;imin(end);imin(end)]);

  y=[y(imin);y(imin(end))+10;y(imin(end))+20];

yoff=mean(y);
xoff=mean(x);
if length(x)>3
F=fit(y-yoff,x-xoff,'poly3');
x=F(y-yoff)+xoff;

elseif length(x)>1
    fitType=['poly',num2str(length(x)-1)];
    F=fit(y-yoff,x-xoff,fitType);
    x=F(y-yoff)+xoff;

end
    


Pbot2=[Pbot;x+50,y];
PbotDist=[0;cumsum(sqrt(sum(diff(Pbot2).^2,2)))];
Ptop2=[Ptop;x+50,y];
PtopDist=[0;cumsum(sqrt(sum(diff(Ptop2).^2,2)))];
%%
%maxDist=max([PtopDist;PbotDist]);
maxDist=400;%fix size of image
interpRange=-100:maxDist+100;
Pbot2=[interp1(PbotDist,Pbot2(:,1),interpRange,'PCHIP')' ...
    interp1(PbotDist,Pbot2(:,2),interpRange,'PCHIP')'];
Pbot2=smooth2a(Pbot2,30,0);
Ptop2=[interp1(PbotDist,Ptop2(:,1),interpRange,'PCHIP')' ...
    interp1(PbotDist,Ptop2(:,2),interpRange,'PCHIP')'];
Ptop2=smooth2a(Ptop2,30,0);

[~,newXbot,newYbot]=wormStraightening((fliplr(Pbot2)),[],150,1);
[~,newXtop,newYtop]=wormStraightening((fliplr(Ptop2)),[],150,1);

%%
midPlaneBuffer=4;
        w=((1:stackSize)-midPlane+midPlaneBuffer)/(midPlaneBuffer*2);
        w=max(w,0);
        w=min(w,1);
        w=reshape(w,1,1,[]);
        [~,~,stackGrid]=meshgrid(1:size(newXtop,2),1:size(newXtop,1),1:stackSize);
        
        newXi=bsxfun(@times,newXbot,(1-w))+bsxfun(@times,newXtop,(w));
        newYi=bsxfun(@times,newYbot,(1-w))+bsxfun(@times,newYtop,(w));
        newXi2=interp3(newXHiFout,newYi,newXi,stackGrid);
        newYi2=interp3(newYHiFout,newYi,newXi,stackGrid);
        segmentStack2=interp3(segmentStack,newYi2,newXi2,stackGrid);
        activityStack2=interp3(activityStack,newYi2,newXi2,stackGrid);

  Iprofile=(smooth(max(smooth2a(nansum(segmentStack2,3),10,10),[],2),30));
  [ymax,imax,ymin,imin]=extrema(Iprofile);
  imin=sort(imin);

  dips=[find(imin<imax(1),1,'last'),find(imin>imax(1),1,'first')];
  dips=imin(dips);
  %%
  if diff(dips)>90 && diff(dips)<140
p=[1:100,linspace(101,dips(1),100),linspace(dips(1)+1,dips(2),110),dips(2)+1:dips(2)+1+290];
segmentStack2=bsxfun(@(A,B) interp1(A,B),segmentStack2,p');
  end

%% plot all images, only works for 2d, probably not here


 
    %%
    
    %    title(num2str(bfIdx));
    metaData(i).metaData.time=time;
    metaData(i).metaData.zVoltage=zPos;
         metaData(i).metaData.iFrame=imageIdx;
         metaData(i).metaData.midPlane=midPlane;
    fileName=['image' num2str(iStack,'%3.5d') '.tif'];
    if saveFlag
        
    tiffwrite([hiResActivityFolder filesep fileName],single(activityStack2),'tif',0);
    tiffwrite([hiResSegmentFolder filesep fileName],single(segmentStack2),'tif',0);
 %   parsave([metaFolder filesep matName],metaData(i),'metaData');
    end
        display(['Finished frame: ' num2str(iStack) ' in ' num2str(toc) 's']);
    catch me
        display(['Error frame: ' num2str(iStack) ' in ' num2str(toc) 's']);
    errors(i).me=me;
    end
    
end
%%
for i=1:length(stackList);
    iStack=stackList(i);
    isempty(metaData(i).metaData)
    if ~isempty(metaData(i).metaData)
            matName=['image' num2str(iStack,'%3.5d')];

    parsave([metaFolder filesep matName],metaData(i),'metaData');
    end
end
