%% track with sameMask

%%
mkdir([dataFolder filesep 'stackDataWhole1Mask'])
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
overwriteFlag=0;
stackList=2:max(stackIdx)-1;

if ~overwriteFlag;
    d=dir([dataFolder filesep 'stackDataWhole1Mask' filesep 'stack*']);
    d={d.name}';
    oldIdx=cell2mat(cellfun(@(x) str2double(x(6:9)),d,'UniformOutput',false));
    stackList=stackList((~ismember(stackList,oldIdx)'));
end
groupSize=100;
%% pick master 
     %% 
     imaster=12;
    imName=['image' num2str(imaster,'%3.5d') '.tif'];
    matName=['image' num2str(imaster,'%3.5d') '.mat'];
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
wormcc=bwconncomp(wormMask);
wormLabelMask=bwlabeln(wormMask,6);
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
%     wormMask=WormSegmentHessian3d_rescale(worm,options);
%     %%
%     %look up intensities on both channels
%     wormLabelMask=bwlabeln(wormMask,6);
stats=regionprops(wormcc,worm,'WeightedCentroid','Area','PixelIdxList');
centroids=reshape([stats.WeightedCentroid],3,[])';

%smooth out activity for interpolating, this is equivalent to expanding the
%ROI's. 
activity=convnfft(activity,gaussianFilter,'same');


Rintensities=cellfun(@(x) mean(worm(worm(x)>quantile(worm(x),.6))),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(worm(x)>quantile(worm(x),.6))),[wormcc.PixelIdxList])';
Rmax=cellfun(@(x) max(worm(x)),[wormcc.PixelIdxList])';
Gmax=cellfun(@(x) max(activity(x)),[wormcc.PixelIdxList])';

    %interpolate Z properly and scale
    realTime=interp1(time,centroids(:,3)); %arb scaling for now
zPlaneIdx=centroids(:,3);
%[~,~,zPix]=cellfun(@(x) ind2sub(size(K),x),{stats.PixelIdxList}','Uniformoutput',false);
%zMax=cellfun(@(x) max(x), zPix);
%zMin=cellfun(@(x) min(x), zPix);

%centroids(:,3)=50*(1+interp1(zPos,centroids(:,3))); %arb scaling for now
Volume=[stats.Area]';

%save outputs in unique file
 outputFile=[dataFolder filesep 'stackDataWhole1Mask' filesep 'stack' num2str(iStack,'%04d') 'data'];
outputData(i).centroids=centroids;
outputData(i).Rintensities=Rintensities;
outputData(i).Gintensities=Gintensities;
outputData(i).Rmax=Rmax;
outputData(i).Gmax=Gmax;
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
        outputFile=[dataFolder filesep 'stackDataWhole1Mask' filesep 'stack' num2str(iStack,'%04d') 'data'];
parsavestruct(outputFile,outputData(i));
        display(['Saving' num2str(iStack,'%04d')]);

    end
end

end


