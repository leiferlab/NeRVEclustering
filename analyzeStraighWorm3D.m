%%  load syncing data
dataFolder='E:\DATA\3DwormData\BrainScanner20140911_182843';
[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);
hiResData=highResTimeTraceAnalysis(dataFolder);


hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
hi2lowTimeOffset=hiResFlashTime-bfFlashTime;
hiResData.frameTime=hiResData.frameTime-hi2lowTimeOffset(1);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

fluor2bfTimeOffset=fluorFlashTime-bfFlashTime;
fluorAll.frameTime=fluorAll.frameTime-fluor2bfTimeOffset(1);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));
Zall=hiResData.Z;
timeAll=hiResData.frameTime;
changes=[0;diff(Zall)<-.6];
zSpikes=[0; diff(changes)<-.5];
stackIdx=cumsum(zSpikes);

lowResFolder=[dataFolder filesep 'lowResFolder'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder'];
metaFolder=[dataFolder filesep  'metaDataFolder'];
%%

options.thresh1=.03; %initial Threshold
options.hthresh=-.0001; %threshold for trace of hessian.
options.minObjSize=500; 
options.maxObjSize=Inf;
options.watershedFilter=1;
options.filterSize=[40,40,40];
options.pad=9;
options.noise=1;
options.show=0;
options.maxSplit=1;
options.minSphericity=.55;
options.valleyRatio=.8;

spikeBuffer=5;
meanSliceFiltLevel=3;
%%
mkdir([dataFolder filesep 'stackDataWhole'])
outputData=struct('FrameData',[]);
outputData=repmat(outputData,max(stackIdx),1);

parfor iStack=7:max(stackIdx)-1
    tic
    try
    imageIdx=find(stackIdx==iStack);
    imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
worm=[];
zPos=[];
time=[];
activity=[];
for iSlice=1:length(imageIdx);
iFrame=imageIdx(iSlice);
        imName=['image' num2str(iFrame,'%3.5d') '.tif'];

        segmentImage=double(imread([hiResSegmentFolder filesep imName],'tif'));
        activityImage=double(imread([hiResActivityFolder filesep imName],'tif'));
        segmentImage=pedistalSubtract(segmentImage);
        activityImage=pedistalSubtract(activityImage);
        worm(:,:,iSlice)=segmentImage;
         activity(:,:,iSlice)=activityImage;
         zPos(iSlice)=Zall(iFrame);
         time(iSlice)=timeAll(iFrame);
end
    
worm(isnan(worm))=0;
goodSlice=squeeze(mean(nanmean(worm)))>meanSliceFiltLevel;
 worm=worm(:,:,goodSlice);
 zPos=zPos(goodSlice);
 time=time(goodSlice);
imageIdx=imageIdx(goodSlice);
%do segmentation
    wormMask=WormSegmentHessian3d_rescale(worm,options);
    
    %look up intensities on both channels
    wormLabelMask=bwlabeln(wormMask,6);
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,'Centroid','Area');
centroids=reshape([stats.Centroid],3,[])';
Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
    realTime=interp1(timeAll,centroids(:,3)); %arb scaling for now

centroids(:,3)=50*(1+interp1(zPos,centroids(:,3))); %arb scaling for now
Volume=[stats.Area]';

%save outputs in unique file
 outputFile=[dataFolder filesep 'stackDataWhole' filesep 'stack' num2str(iStack,'%04d') 'data'];
outputData(iStack).centroids=centroids;
outputData(iStack).Rintensities=Rintensities;
outputData(iStack).Gintensities=Gintensities;
outputData(iStack).Volume=Volume;
outputData(iStack).time=realTime;
outputData(iStack).wormMask=wormLabelMask;

outputData(iStack).FrameData.imageIdx=imageIdx;
outputData(iStack).FrameData.zPos=zPos;
 outputData(iStack).FrameData.time=time;


parsavestruct(outputFile,outputData(iStack));
display(['Completed stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);

    catch
    display(['Error in stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);
    
    end
    


    
    
end


%%


