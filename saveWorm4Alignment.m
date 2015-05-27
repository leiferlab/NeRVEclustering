alignFolder=[dataFolder filesep 'wormAlignmentImages'];
mkdir(alignFolder)


[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);


rows=1200;
cols=600;
nPix=rows*cols;


bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


%% set up low magvideos

aviFiles=dir([dataFolder filesep '20*.avi']);
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

%%

for iStack=hasPoints(1:50)';
iStack
try
    %%
Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);
behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);



   segmentChannel2=[]; segmentChannel3=[]; segmentChannel4=[]; segmentChannel5=[];

    hiResIdx=find(hiResData.stackIdx==iStack)+ zOffset;
    zRange=hiResData.Z(hiResIdx-zOffset);
    
    
          bfIdx=round(bfIdxLookup(hiResIdx));
            fluorIdx=round(fluorIdxLookup(hiResIdx));
         fluorIdxRange=[min(fluorIdx) max(fluorIdx)];
        bfIdxRange=[min(bfIdx) max(bfIdx)];
[~,ia,ib]=unique(fluorIdx);
         fluorFrame=read(fluorVidObj,fluorIdxRange);
            bfFrame = read(behaviorVidObj,bfIdxRange);
            fluorFrame=squeeze(fluorFrame(:,:,1,:));
%             bfFrame=squeeze(bfFrame(:,:,1,:));
%                   bfFrame=imwarp(bfFrame,invert(lowResFluor2BF.t_concord),...
%                 'OutputView',lowResFluor2BF.Rsegment);
%               bfFrame=imwarp(bfFrame,Hi2LowResF.t_concord,...
%                 'OutputView',Hi2LowResF.Rsegment);
%                bfFrame=bfFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
% 
%             fluorFrame2=imwarp(fluorFrame,Hi2LowResF.t_concord,...
%                 'OutputView',Hi2LowResF.Rsegment);
%             fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
            
 fluorFrame=double(fluorFrame(:,:,ib));           
 %bfFrame=double(bfFrame(:,:,ib));           
            
    status=fseek(Fid,2*(hiResIdx(1))*nPix,-1);
    pixelValues=fread(Fid,nPix*(length(hiResIdx)),'uint16',0,'l');
 
    
    hiResImage=reshape(pixelValues,rows,cols,length(hiResIdx));
    
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
    activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
    segmentChannel=pedistalSubtract(segmentChannel);
    activityChannel=pedistalSubtract(activityChannel);
    
    alignImageHi=[alignFolder filesep 'hiRes' num2str(iStack,'%3.5d')];
    alignImageLo=[alignFolder filesep 'lowRes' num2str(iStack,'%3.5d')];
    
    tiffwrite(alignImageHi,max(segmentChannel(:,:,18:22),[],3),'.tif');
    tiffwrite(alignImageLo,mean(fluorFrame(:,:,18:22),3),'.tif');
    display(['Finished '  num2str(iStack,'%3.5d')])
catch ME
    ME
end

end
