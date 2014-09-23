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



%% load centerline data
centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile=centerLineFile.name;
centerline=load([dataFolder filesep centerLineFile],'centerline');
centerline=centerline.centerline;

%% make lookup tables for indices
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'PCHIP');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');

%%
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder'];
hiResSegmentFolder=[dataFolder filesep 'hiResSegmentFolder'];

imageFiles=dir([hiResActivityFolder filesep '*.tif']);
imageFiles={imageFiles.name}';


%% segment subimages and create masks
paraPool=gcp('nocreate');
if ~isempty(paraPool)
paraPool=parpool('local',8,'IdleTimeout', 600);
else
    delete(paraPool)
    paraPool=parpool('local',8,'IdleTimeout', 600);

end
    
%%
outputData=struct('imName',[]);
outputData=repmat(outputData,length(imageFiles),1);
parfor iImage=1:length(imageFiles);
    iFrame=str2double(imageFiles{iImage}(6:11)); %pull out frame idx, must be in for image*****.tif
    worm=imread([hiResSegmentFolder filesep imageFiles{iImage}],'tif');
    activity=imread([hiResActivityFolder filesep imageFiles{iImage}],'tif');

  wormMask=WormSegmentHessian2D_whole(worm);
   wormMask= bwmorph(wormMask,'clean');
   
       %look up intensities on both channels
wormLabelMask=bwlabeln(wormMask);
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,'Centroid','Area');
centroids=cell2mat({stats.Centroid}');
Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';
Volume=[stats.Area]';


outputFile=[dataFolder filesep 'stackData2D' filesep 'stack' num2str(iFrame,'%05d') 'data'];
 outputData(iImage).centroids=centroids;
 outputData(iImage).Rintensities=Rintensities;
outputData(iImage).Gintensities=Gintensities;
 outputData(iImage).Volume=Volume;
 outputData(iImage).time=hiResData.frameTime(iFrame);
 outputData(iImage).wormMask=wormLabelMask;
 outputData(iImage).imName=imageFiles{iImage};
 
outputData(iImage)
% parsavestruct(outputFile,outputData(iImage));
end



%%
parfor iImage=1:length(imageFiles);
tic
    iFrame=str2double(imageFiles{iImage}(6:11)); %pull out frame idx, must be in for image*****.tif
 if iFrame<12865    % only going up to when Z scanning starts
     
 
 
 %progressbar((iFrame-1750)/12865)
 
    worm=imread([hiResSegmentFolder filesep imageFiles{iImage}],'tif');
    activity=imread([hiResActivityFolder filesep imageFiles{iImage}],'tif');

  wormMask=WormSegmentHessian2D_whole(worm);
   wormMask= bwmorph(wormMask,'clean');
   
       %look up intensities on both channels
wormLabelMask=bwlabeln(wormMask);
wormcc=bwconncomp(wormMask);
stats=regionprops(wormcc,'Centroid','Area');
centroids=cell2mat({stats.Centroid}');
Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';
Volume=[stats.Area]';

outputFile=[dataFolder filesep 'stackData2D' filesep 'stack' num2str(iFrame,'%05d') 'data'];
outputData(iImage).centroids=centroids;
outputData(iImage).Rintensities=Rintensities;
outputData(iImage).Gintensities=Gintensities;
outputData(iImage).Volume=Volume;
outputData(iImage).time=hiResData.frameTime(iFrame);
outputData(iImage).wormMask=wormLabelMask;
outputData(iImage).imName=imageFiles{iImage};
%outputData(iImage)
parsavestruct(outputFile,outputData(iImage));

   display(['Analysis for frame ' imageFiles{iImage} ' completed in ' num2str(toc) ' seconds']);
   
 end
end


   



