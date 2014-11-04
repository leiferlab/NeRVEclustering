%%  load syncing data
dataFolder=uipickfiles;
dataFolder=dataFolder{1};
[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));
end

hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
[~,Hi2bf]=flashTimeAlign2(hiResFlashTime,bfFlashTime);
flashDiff=hiResFlashTime(Hi2bf)-bfFlashTime;
flashDiff=flashDiff-min(flashDiff);
f_hiResTime=fit(hiResFlashTime(Hi2bf),bfFlashTime,'poly1','Weight',exp(-flashDiff.^2));

hiResData.frameTime=f_hiResTime(hiResData.frameTime);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

[~,bf2fluor]=flashTimeAlign2(bfFlashTime,fluorFlashTime);

flashDiff=fluorFlashTime-bfFlashTime(bf2fluor);
flashDiff=flashDiff-min(flashDiff);

f_fluorTime=fit(fluorFlashTime,bfFlashTime(bf2fluor),'poly1','Weight',exp(-flashDiff.^2));

if f_fluorTime.p1<.1
    f_fluorTime.p1=1;
end


fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));

%%
firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));
Zall=hiResData.Z;
timeAll=hiResData.frameTime;
stackIdx=hiResData.stackIdx;

lowResFolder=[dataFolder filesep 'lowResFolder'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder3D'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder3D'];
metaFolder=[dataFolder filesep  'metaDataFolder3D'];
%%

options.thresh1=.03; %initial Threshold
options.hthresh=-.001; %threshold for trace of hessian.
options.minObjSize=140; 
options.maxObjSize=Inf;
options.watershedFilter=1;
options.filterSize=[5 5 3];
options.pad=9;
options.noise=1;
options.show=0;
options.maxSplit=1;
options.minSphericity=.55;
options.valleyRatio=.8;
options.scaleFactor=[1,1,3];
spikeBuffer=0;
meanSliceFiltLevel=3;

gaussianFilter=fspecial('gaussian',[30,30],5);
gaussianFilter=convnfft(gaussianFilter,permute(gausswin(6,2),[2,3,1]));
gaussianFilter=gaussianFilter/sum(gaussianFilter(:));


%%
mkdir([dataFolder filesep 'stackDataWhole'])
outputData=[];
outputData.centroids=[];
outputData.Rintensities=[];
outputData.Gintensities=[];
outputData.Volume=[];
outputData.time=[];
outputData.wormMask=[];
outputDat.zPlaneIdx=[];
outputData.zMax=[];
outputData.zMin=[];
output0=outputData;
overwriteFlag=1;
stackList=2:max(stackIdx)-1;

if ~overwriteFlag;
    d=dir([dataFolder filesep 'stackDataWhole' filesep 'stack*']);
    d={d.name}';
    oldIdx=cell2mat(cellfun(@(x) str2double(x(6:9)),d,'UniformOutput',false));
    stackList=stackList((~ismember(stackList,oldIdx)'));
end
groupSize=100;
%%
for iGroup=1:length(stackList)/groupSize

group=groupSize*(iGroup-1):groupSize*(iGroup);

group=group(group<length(stackList));
group=group(group>1);

subStackList=stackList(group);
outputData=repmat(output0,length(subStackList),1);

parfor i=1:length(group);

    try
        
        iStack=subStackList(i);
    tic    
        %%
    imName=['image' num2str(iStack,'%3.5d') '.tif'];
    matName=['image' num2str(iStack,'%3.5d') '.mat'];
    imageIdx=find(stackIdx==iStack);
    imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
    metaData=load([metaFolder filesep matName]);
    metaData=metaData.metaData.metaData;
    
worm=stackLoad([hiResSegmentFolder filesep imName]);
worm=pedistalSubtract(worm);
worm(isnan(worm))=0;
zPos=metaData.zVoltage;
time=metaData.time;
activity=stackLoad([hiResActivityFolder filesep imName]);
activity=pedistalSubtract(activity);
activity(isnan(activity))=0;

%% do segmentation
    wormMask=WormSegmentHessian3d_rescale(worm,options);
    %%
    %look up intensities on both channels
    wormLabelMask=bwlabeln(wormMask,6);
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,worm,'WeightedCentroid','Area','PixelIdxList');
centroids=reshape([stats.WeightedCentroid],3,[])';

%smooth out activity for interpolating, this is equivalent to expanding the
%ROI's. 
activity=convnfft(activity,gaussianFilter,'same');


Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
    realTime=interp1(time,centroids(:,3)); %arb scaling for now
zPlaneIdx=centroids(:,3);
%[~,~,zPix]=cellfun(@(x) ind2sub(size(K),x),{stats.PixelIdxList}','Uniformoutput',false);
%zMax=cellfun(@(x) max(x), zPix);
%zMin=cellfun(@(x) min(x), zPix);

centroids(:,3)=50*(1+interp1(zPos,centroids(:,3))); %arb scaling for now
Volume=[stats.Area]';

%save outputs in unique file
 outputFile=[dataFolder filesep 'stackDataWhole' filesep 'stack' num2str(iStack,'%04d') 'data'];
outputData(i).centroids=centroids;
outputData(i).Rintensities=Rintensities;
outputData(i).Gintensities=Gintensities;
outputData(i).Volume=Volume;
outputData(i).time=time;
ouptutData(i).centroidTime=realTime;
outputData(i).wormMask=wormLabelMask;
outputData(i).zPlaneIdx=zPlaneIdx;
%outputData(i).zMax=zMax;
%outputData(i).zMin=zMin;
%parsavestruct(outputFile,outputData(i));
display(['Completed stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);

    catch ME
    display(['Error in stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);
    ME
    end
    
end
%%
for i=1:length(outputData);
    if ~isempty(outputData(i).centroids)
        iStack=subStackList(i);
        outputFile=[dataFolder filesep 'stackDataWhole' filesep 'stack' num2str(iStack,'%04d') 'data'];
parsavestruct(outputFile,outputData(i));
        display(['Saving' num2str(iStack,'%04d')]);

    end
end

end



%%


