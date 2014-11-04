Ztype='traingle';
imSize=[1200,600];
row=imSize(1);
col=imSize(2);
nPix=row*col;
%%  load syncing data
dataFolder='E:\DATA\3DwormData\20141020\20141020_1610Data';


[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);
if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
hiResData=highResTimeTraceAnalysisTriangle3(dataFolder,imSize(1),imSize(2));
end

hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
hi2lowTimeOffset=hiResFlashTime(1)-bfFlashTime(1);
hiResData.frameTime=hiResData.frameTime-hi2lowTimeOffset(1);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

fluor2bfTimeOffset=fluorFlashTime(1)-bfFlashTime(5);
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
centerline=load([dataFolder filesep centerLineFile{1}],'centerline');
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
if movieFlag
    plotFlag=1;
vidOut=VideoWriter([dataFolder filesep 'HiMagOnly.avi']);
vidOut.FrameRate=20;

open(vidOut);
end
spikeBuffer=5;
meanSliceFiltLevel=4;
yWindow=300;
xWindow=200;
hiResCorrection=[0,0];
splitFlag=0;
stackIdx=hiResData.stackIdx;
%%
%progressbar(0,0)

testAlignmentFolder=[dataFolder filesep 'testAlignmentFolder' filesep];
mkdir(testAlignmentFolder);
for iStack=2:max(stackIdx);
    
        imageIdx=find(stackIdx==iStack);
    imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
zPos=hiResData.Z(imageIdx);
time=hiResData.frameTime(imageIdx);
iSlice=(round(length(imageIdx)/2));
        iFrame=imageIdx(iSlice);

    %%
      %  progressbar(iStack/max(stackIdx),iSlice/length(imageIdx))
    iTime=hiResData.frameTime(iFrame);
    %interpolate using time to get low res idx 
    bfIdx=bfIdxLookup(iFrame);
fluorIdx=fluorIdxLookup(iFrame);
hiResIdx=(iFrame);
  status=fseek(Fid,2*hiResIdx*nPix,-1);

  pixelValues=fread(Fid,nPix,'uint16',0,'l');

hiResImage=(reshape(pixelValues,row,col));
hiResImage=hiResImage-backgroundImage;
hiResImage(hiResImage<0)=0;
 activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
 segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
 activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
activityChannel=pedistalSubtract(activityChannel,5);
% segmentStack(:,:,iSlice)=segmentChannel;
% activityStack(:,:,iSlice)=activityChannel;
% % 
%  hiResImage2=imwarp(hiResImage,Hi2LowRes.t_concord,...
%     'OutputView',Hi2LowRes.Rsegment);

fluorFrame=read(fluorVidObj,round(fluorIdx));
bfFrame = read(behaviorVidObj,round(bfIdx));
fluorFrame=fluorFrame(:,:,1);
bfFrame=bfFrame(:,:,1);

 fluorFrame=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
     'OutputView',lowResFluor2BF.Rsegment);
fluorFrame=imwarp(fluorFrame,Hi2LowRes.t_concord,  'OutputView',Hi2LowRes.Rsegment);
 fluorFrame=fluorFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

  bfFrame=imwarp(bfFrame,Hi2LowRes.t_concord,...
    'OutputView',Hi2LowRes.Rsegment);
 bfFrame=bfFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

subplot(1,3,1);
imagesc(fluorFrame);
axis equal
subplot(1,3,2);
imagesc(segmentChannel);
axis equal
subplot(1,3,3)
imagesc(bfFrame);
axis equal
drawnow;

% tiffwrite([testAlignmentFolder 'hiRes' num2str(iStack)],hiResImage,'tif',0);
% tiffwrite([testAlignmentFolder 'fluorFrame' num2str(iStack)],fluorFrame,'tif',0);

end
