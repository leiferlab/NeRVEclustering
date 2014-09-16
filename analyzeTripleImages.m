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

%works amazingly well, flash frames off by less than .05 seconds. 

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
    movies=uipickfiles;
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

behaviorVidObj = VideoReader(behaviorMovie);
 fluorVidObj= VideoReader(fluorMovie);

Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'nearest');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'nearest');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));

%%
for iFrame=firstFullFrame:length(hiResData.frameTime);
    iTime=hiResData.frameTime(iFrame);
    %interpolate using time to get low res idx 
    bfIdx=bfIdxLookup(iFrame);
fluorIdx=fluorIdxLookup(iFrame);
hiResIdx=hiResLookup(iFrame);
  status=fseek(Fid,2*hiResIdx*1024^2,-1);

  pixelValues=fread(Fid,1024^2,'uint16',0,'l');

hiResImage=rot90(reshape(pixelValues,1024,1024)',2);

fluorFrame=read(fluorVidObj,fluorIdx);
bfFrame = read(behaviorVidObj,bfIdx);
    
    subplot(2,2,[1,3]);
    imagesc(hiResImage);
    title(num2str(hiResIdx));
    subplot(2,2,2);
    imagesc(fluorFrame);
    title(num2str(fluorIdx));
    subplot(2,2,4);
    imagesc(bfFrame);
        title(num2str(bfIdx));

    drawnow;
    
end
