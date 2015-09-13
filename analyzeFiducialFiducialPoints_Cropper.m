%%  load syncing data

dataFolder=uipickfiles;
dataFolder=dataFolder{1};
%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
%%
imSize=[1200 600];
nPix=prod(imSize);
if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));
end

%z2ImageIdxOffset=-9;


%%

% bfIdxList=1:length(bfAll.frameTime);
% fluorIdxList=1:length(fluorAll.frameTime);
% bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
% fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');
% 
% [hiImageIdx,ib]=unique(hiResData.imageIdx);
% hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));
% 
% %%
% firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
% firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));
Zall=hiResData.Z;
timeAll=hiResData.frameTime;
stackIdx=hiResData.stackIdx;

lowResFolder=[dataFolder filesep 'lowResFoldertest'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder3Dtest'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder3Dtest'];
metaFolder=[dataFolder filesep  'metaDataFolder3Dtest'];

%% load Fiducials file
fiducialFile=dir([dataFolder filesep '*iducial*']);
fiducialFile={fiducialFile.name}';
if length(fiducialFile)~=1
        display('Select model file');

    fiducialFile=uipickfiles('FilterSpec',dataFolder);
    fiducialFile=load(fiducialFile{1});
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;
else
    fiducialFile=load([dataFolder filesep fiducialFile{1}]);
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;

end



%%
try
    alignments=load([dataFolder filesep 'alignments.mat']);
    alignments=alignments.alignments;
    S2AHiRes=alignments.S2AHiRes;
    rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
catch
display('Select Hi Res Alignment')

S2AHiResFile=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiResFile{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;
end

%%
gaussianFilter=fspecial('gaussian',[5,5],2);
gaussianFilter=convnfft(gaussianFilter,permute(gausswin(6,1),[2,3,1]));
gaussianFilter=gaussianFilter/sum(gaussianFilter(:));


%%
%mkdir([dataFolder filesep 'stackDataFiducials'])
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

R=10; %radius of dilation around each point
if ~overwriteFlag;
    d=dir([dataFolder filesep 'stackDataWhole' filesep 'stack*']);
    d={d.name}';
    oldIdx=cell2mat(cellfun(@(x) str2double(x(6:9)),d,'UniformOutput',false));
    stackList=stackList((~ismember(stackList,oldIdx)'));
end
groupSize=100;

%% make label mask based on reference points and image size. 

% Fpoints=sub2ind(size(refstack),round(masterFiducials(:,1)),...
%     round(masterFiducials(:,2)),round(masterFiducials(:,3)));
% 
% wormLabelMask=zeros(size(refstack));
% wormLabelMask(Fpoints)=1:length(Fpoints);
% Temp=fspecial('gaussian',[25 25],8);
% Temp=convnfft(Temp,permute(gausswin(6,2),[2,3,1]));
% Temp=Temp/sum(Temp(:));
% 
% wormLabelMask=imdilate(wormLabelMask,Temp>max(Temp(:)/3));
% wormcc=bwconncomp(wormLabelMask);
% stats=regionprops(wormLabelMask,'Centroid','Area');
% centroids=reshape([stats.Centroid],3,[])';




if 0
    
      
    greenBackground=gback((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    redBackground=rback((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
    redBackground=imwarp(redBackground,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
greenBackground=greenBackground-90;
greenBackground=greenBackground/mean(greenBackground(:));
greenBackground(greenBackground<.1)=.1;
redBackground=redBackground-90;
redBackground=redBackground/mean(redBackground(:));
redBackground(redBackground<.1)=.1;

else
    greenBackground=1;
    redBackground=1;
    
    
end

%%


annotated=find(cell2mat(cellfun(@(x) ~isempty(cell2mat(x(1))),fiducialPoints,'uniformOutput',0)));
%%
progressbar(0);
boxR=[10,10,1];
[boxX,boxY,boxZ]=meshgrid(-boxR(1):boxR(1),-boxR(2):boxR(2),-boxR(3):boxR(3));
Rcrop=zeros([size( boxX) length(annotated) size(fiducialPoints{1},1)]);
Gcrop=Rcrop;

for iStack=1:length(annotated)
  %  progressbar(iStack/length(annotated));

      currentFiducials=fiducialPoints{annotated(iStack)};
       Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

        fiducialIdx=find(cell2mat((cellfun(@(x)( ~isempty(x)),currentFiducials(:,1),'uniformoutput',0))));
        fiducialIdx(isnan(cell2mat(currentFiducials(fiducialIdx,1))))=[];
          fiducialZ=round(cell2mat(currentFiducials(fiducialIdx,4)));
    %    [fiducialPlanes,ia,ib]=unique(fiducialZ);

        currentFiducials=cell2mat(currentFiducials(fiducialIdx,1:4));    
        Rout=zeros(12,1);
        Gout=Rout;
        backgroundOut=Gout;
        
        rFiducialPoints=round(currentFiducials(:,[1,2,4]));
        
%%

lowLimit=min(rFiducialPoints(:,3)-boxR(3));
upperLimit=max(rFiducialPoints(:,3)+boxR(3));
nSlices=upperLimit-lowLimit+1;
   %     hiResIdx=(fiducialZ(i)-boxR(3));

%    inPlaneFiducials=currentFiducials(ib==i,:);
    
    status=fseek(Fid,2*(lowLimit)*nPix,-1);
    
    pixelValues=fread(Fid,nPix*nSlices,'uint16',0,'l');
    

    hiResImage=(reshape(pixelValues,imSize(1),imSize(2),nSlices));
    
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
    activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
%%
for i=1:size(rFiducialPoints,1)
    %%
progressbar(iStack/length(annotated),i/size(rFiducialPoints,1));
%     segmentChannel=bsxfun(@rdivide,segmentChannel,greenBackground);
%     activityChannel=bsxfun(@rdivide,activityChannel,redBackground);
%  
currentPoint=rFiducialPoints(i,:);
currentPoint(3)=currentPoint(3)-lowLimit+boxR(3);
zRange=currentPoint(3)+(-boxR(3):boxR(3));
    [yspace,xspace]=meshgrid(currentPoint(2)+(-boxR(1):boxR(1))...
    ,currentPoint(1)+(-boxR(2):boxR(2)));
xspacemap=xspace<(rect1(3)-rect1(1))& xspace>1;
yspacemap=yspace<(rect1(4)-rect1(2)) & yspace>1;
inImageMap=and(xspacemap, yspacemap);
inImageMap2=repmat(inImageMap',[1,1,boxR(3)*2+1]);
Rtemp=nan(2*boxR+1);
Gtemp=Rtemp;

Rtemp(inImageMap2)=segmentChannel(yspace(1, any(inImageMap,1)),xspace(any(inImageMap,2),1),zRange);
Gtemp(inImageMap2)=activityChannel(yspace(1, any(inImageMap,1)),xspace(any(inImageMap,2),1),zRange);

% Rtemp=segmentChannel(rFiducialPoints(i,2)+(-boxR(1):boxR(1))...
%     ,rFiducialPoints(i,1)+(-boxR(2):boxR(2)),:);
% Gtemp=activityChannel(rFiducialPoints(i,2)+(-boxR(1):boxR(1))...
%     ,rFiducialPoints(i,1)+(-boxR(2):boxR(2)),:);
%RBW=(normalizeRange(Rtemp)>graythresh(normalizeRange(Rtemp)));

Rcrop(:,:,:,iStack,fiducialIdx(i))=Rtemp;
Gcrop(:,:,:,iStack,fiducialIdx(i))=Gtemp;

end
fclose(Fid);
% Rfiducials{iStack}=Rout;
% Gfiducials{iStack}=Gout;
% BackgroundF{iStack}=backgroundOut;
end
%%

SliceBrowser((squeeze(mean(mean(Gcrop(:,:,:,:,:),2),1))));

SliceBrowser((squeeze(mean(mean(Rcrop(:,:,:,:,:),3),1))));
%%
Gcat=[];
Rcat=[];
Gcatxy=[];Gproj=[];
Rcatxy=[];Rproj=[];

for i=1:size(Rcrop,5)
    
    Gtemp=(smooth2a(squeeze(mean(mean(Gcrop(:,:,:,:,i),2),1)),2,1));
    Rtemp=(smooth2a(squeeze(mean(mean(Rcrop(:,:,:,:,i),2),1)),2,1));
    Gtempxy=(smooth2a(squeeze(mean(mean(Gcrop(:,:,:,:,i),3),2)),2,1));
    Rtempxy=(smooth2a(squeeze(mean(mean(Rcrop(:,:,:,:,i),3),2)),2,1));
        
    
    Rtemp(Rtemp>500)=nan;
Gtemp(Gtemp>500)=nan;
Gcat=cat(1,Gcat,normalizeRange(Gtemp(:,10:end)));
Rcat=cat(1,Rcat,normalizeRange(Rtemp(:,10:end)));
Gcatxy=cat(1,Gcatxy,normalizeRange(Gtempxy(:,10:end)));
Rcatxy=cat(1,Rcatxy,normalizeRange(Rtempxy(:,10:end)));
Rproj=cat(1,Rproj,mean(Rtemp(:,10:end)));
Gproj=cat(1,Gproj,mean(Gtemp(:,10:end)));
end



%%


%%
save([dataFolder filesep 'wormFiducialIntensities_bot'],'Rcrop','Gcrop');

