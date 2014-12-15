%%  load syncing data
imSize=[1200 600];
dataFolder=uipickfiles;
dataFolder=dataFolder{1};
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
z2ImageIdxOffset=-9;


%% load centerline data
refStackFile=dir([dataFolder filesep '*refStack*']);
refStackFile={refStackFile.name}';
if length(refStackFile)~=1
    refStackFile=uipickfiles('FilterSpec',dataFolder);
    load(refStackFile{1});
else
    load([dataFolder filesep refStackFile{1}]);
end
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

lowResFolder=[dataFolder filesep 'lowResFoldertest'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder3Dtest'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder3Dtest'];
metaFolder=[dataFolder filesep  'metaDataFolder3Dtest'];
%%

options.thresh1=.03; %initial Threshold
options.hthresh=-.001; %threshold for trace of hessian.
options.minObjSize=140; 
options.maxObjSize=Inf;
options.watershedFilter=1;
options.filterSize=[15 15 3];
options.pad=9;
options.noise=1;
options.show=0;
options.maxSplit=1;
options.minSphericity=.55;
options.valleyRatio=.8;
options.scaleFactor=[1,1,6];
spikeBuffer=0;
meanSliceFiltLevel=3;

gaussianFilter=fspecial('gaussian',[5,5],2);
gaussianFilter=convnfft(gaussianFilter,permute(gausswin(6,1),[2,3,1]));
gaussianFilter=gaussianFilter/sum(gaussianFilter(:));


%%
mkdir([dataFolder filesep 'stackDataFiducials'])
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

%% make label mask based on reference points and image size. 

Fpoints=sub2ind(size(refstack),round(masterFiducials(:,1)),...
    round(masterFiducials(:,2)),round(masterFiducials(:,3)));

wormLabelMask=zeros(size(refstack));
wormLabelMask(Fpoints)=1:length(Fpoints);
Temp=fspecial('gaussian',[25 25],8);
Temp=convnfft(Temp,permute(gausswin(6,2),[2,3,1]));
Temp=Temp/sum(Temp(:));

wormLabelMask=imdilate(wormLabelMask,Temp>max(Temp(:)/3));
wormcc=bwconncomp(wormLabelMask);
stats=regionprops(wormLabelMask,'Centroid','Area');
centroids=reshape([stats.Centroid],3,[])';

%%
outputAll=[];
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
%worm=pedistalSubtract(worm);
worm(isnan(worm))=0;
zPos=metaData.zVoltage;
time=metaData.time;
activity=stackLoad([hiResActivityFolder filesep imName]);
%activity=pedistalSubtract(activity);
activity(isnan(activity))=0;

%% do segmentation
    %wormMask=WormSegmentHessian3d_rescale(worm,options);
    %%
    %look up intensities on both channels
   % wormLabelMask=bwlabeln(wormMask,6);

%smooth out activity for interpolating, this is equivalent to expanding the
%ROI's. 
activity=convnfft(activity,gaussianFilter,'same');
Rall=accumarray(wormLabelMask(wormLabelMask>0),worm(wormLabelMask>0),[],...
    @(x) {x});
Gall=accumarray(wormLabelMask(wormLabelMask>0),activity(wormLabelMask>0),[],...
    @(x) {x});

Rintensities=cellfun(@(x) trimmean(x,20), Rall,'uniformoutput',0);
Gintensities=cellfun(@(a,b) mean(a(b>max(b(:)/2)),Gall,Rall,'uniformoutput',0));
Rintensities=cell2mat(Rintensities);Gintensities=cell2mat(Gintensities);
% Rmax=cellfun(@(x) max(worm(x)),[wormcc.PixelIdxList])';
% Gmax=cellfun(@(x) max(activity(x)),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
    realTime=interp1(time,centroids(:,3)); %arb scaling for now
zPlaneIdx=centroids(:,3);
%[~,~,zPix]=cellfun(@(x) ind2sub(size(K),x),{stats.PixelIdxList}','Uniformoutput',false);
%zMax=cellfun(@(x) max(x), zPix);
%zMin=cellfun(@(x) min(x), zPix);

%centroids(:,3)=50*(1+interp1(zPos,centroids(:,3))); %arb scaling for now
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
        outputFile=[dataFolder filesep 'stackDataFiducials' filesep 'stack' num2str(iStack,'%04d') 'data'];
parsavestruct(outputFile,outputData(i));
        display(['Saving' num2str(iStack,'%04d')]);

    end
end
outputAll=[outputAll;outputData];

end





%%


