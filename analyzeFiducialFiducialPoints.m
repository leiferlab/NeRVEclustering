%%  load syncing data
imSize=[1200 600];
nPix=prod(imSize);
dataFolder=uipickfiles;
dataFolder=dataFolder{1};
%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));
end

z2ImageIdxOffset=-9;


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
else
    fiducialFile=load([dataFolder filesep fiducialFile{1}]);
    fiducialPoints=fiducialFile.fiducialPoints;
end



%%
display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;


%%
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

%%

annotated=find(cell2mat(cellfun(@(x) ~isempty(cell2mat(x(1))),fiducialPoints,'uniformOutput',0)));
%%
progressbar(0,0);
for iStack=annotated'
      currentFiducials=fiducialPoints{iStack};
        Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

        fiducialIdx=find(cell2mat((cellfun(@(x)( ~isempty(x)),currentFiducials(:,1),'uniformoutput',0))));
        fiducialIdx(isnan(cell2mat(currentFiducials(fiducialIdx,1))))=[];
          fiducialZ=round(cell2mat(currentFiducials(fiducialIdx,4)));
       fiducialZ=fiducialZ;
        [fiducialPlanes,ia,ib]=unique(fiducialZ);

        currentFiducials=cell2mat(currentFiducials(fiducialIdx,:));    
        Rout=zeros(12,1);
        Gout=Rout;
        backgroundOut=Gout;
%%

for i=1:length(fiducialPlanes);
progressbar(iStack/max(annotated),i/length(fiducialPlanes));

        hiResIdx=(fiducialPlanes(i));

    inPlaneFiducials=currentFiducials(ib==i,:);
    
    status=fseek(Fid,2*(hiResIdx)*nPix,-1);
    
    pixelValues=fread(Fid,nPix,'uint16',0,'l');
    

    hiResImage=(reshape(pixelValues,imSize(1),imSize(2)));
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
    activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
    activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
    
    Fpoints=sub2ind(size(segmentChannel),round(inPlaneFiducials(:,2)),...
    round(inPlaneFiducials(:,1)));

wormLabelMask=zeros(size(segmentChannel));
wormLabelMask(Fpoints)=1:length(Fpoints);
wormLabelMask=imdilate(wormLabelMask, strel('disk', R));
% figure(2)
% imagesc(wormLabelMask)
% hold on
% contour(wormLabelMask,[0,0]);
% scatter(inPlaneFiducials(:,1),inPlaneFiducials(:,2),'rx');
% hold off
% pause(1)

Rall=accumarray(wormLabelMask(wormLabelMask>0),segmentChannel(wormLabelMask>0),[],...
    @(x) {x});
Gall=accumarray(wormLabelMask(wormLabelMask>0),activityChannel(wormLabelMask>0),[],...
    @(x) {x});

Rintensities=cellfun(@(x) trimmean(x,20), Rall,'uniformoutput',0);
Gintensities=cellfun(@(a,b) mean(a(b>max(b(:)/2))),Gall,Rall,'uniformoutput',0);
Rintensities=cell2mat(Rintensities);Gintensities=cell2mat(Gintensities);

Rout( fiducialIdx(ib==i))=Rintensities;
Gout(fiducialIdx(ib==i))=Gintensities;
backgroundOut(fiducialIdx(ib==i))=mean(segmentChannel(segmentChannel>quantile(segmentChannel(:),.8)));
% Rmax=cellfun(@(x) max(worm(x)),[wormcc.PixelIdxList])';
% Gmax=cellfun(@(x) max(activity(x)),[wormcc.PixelIdxList])';

end
fclose(Fid)
Rfiducials{iStack}=Rout;
Gfiducials{iStack}=Gout;
BackgroundF{iStack}=backgroundOut;
end
%%

Gall=cell2mat(Gfiducials);
Rall=cell2mat(Rfiducials);
BackgroundAll=cell2mat(BackgroundF);
Rall(Rall==0)=nan;
Gall(Gall==0)=nan;
Rall=Rall(~all(isnan(Rall),2),:);
Gall=Gall(~all(isnan(Gall),2),:);
%%
for i=1:size(Rall,1)
    figure
    plot((Rall(i,:)),'r')
hold on
plot((Gall(i,:)),'g')
time=(1:size(Rall,2))';
I=Rall(i,:)';
time=time(~isnan(I));
I=I(~isnan(I));

f=fit(time,I,'poly1');
plot(f)



 figure
 plot((Gall(i,:)./nanmean(Gall(i,:)))./(Rall(i,:)./f(1:size(Rall,2))'));
% plot((Rall(i,:)./f(1:size(Rall,2))'));

%plot(2*BackgroundAll(i,:),'black')

end

%%


%%
save([dataFolder filesep 'wormFiducialIntensities'],'Rall','Gall');

