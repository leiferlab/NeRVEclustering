display('Select submission Folder')
submissionFolder=uipickfiles('filterSpec','Y:\Jeff\');

submissionFolder=submissionFolder{1};
fileListAll=dir([submissionFolder filesep 'bot*.mat']);
display('Select pointstats')
pointStatsFile=uipickfiles('filterspec',fileparts(submissionFolder));

%%
load(pointStatsFile{1});
pointStatsNew=pointStats2;

%% select image folder

imageFolder=uipickfiles('filterSpec','V:');
imageFolder=imageFolder{1};
dataFolder=fileparts(imageFolder)

%% 
    iPoint=str2double(fileListAll(end).name(11:15));
    data=load([submissionFolder filesep fileListAll(10).name]);
    
    xmean=nanmean(data.comparePointEstimate_x)';
    
  newXAll=nan(iPoint,length(xmean));
  newYAll=newXAll;
  newZAll=newXAll;
  oldXAll=newXAll;
  oldYAll=newXAll;
  oldZAll=newXAll;
  zScoreAll=newXAll;
  compareAllX=zeros(size(data.comparePointEstimate_x,1),...
      size(data.comparePointEstimate_x,2),length(fileListAll)/3);
  compareAllY=compareAllX;
  compareAllZ=compareAllX;
  
for i=1:length(fileListAll)
   % i=pointStats2(iPS).stackIdx;
    fileName=['botChecker' num2str(i,'%3.5d') '*'];
    display(['loading ' fileName]);
    fileList=dir([submissionFolder filesep fileName]);
    fileList={fileList.name};
 clear data
 if ~isempty(fileList)
for j=1:length(fileList)
    data(j)=load([submissionFolder filesep fileList{j}]);
end
comparePointEstimate_x=cell2mat(permute({data.comparePointEstimate_x},[1 3 2]));
comparePointEstimate_x=nansum(comparePointEstimate_x,3);
comparePointEstimate_x(comparePointEstimate_x==0)=nan;

comparePointEstimate_y=cell2mat(permute({data.comparePointEstimate_y},[1 3 2]));
comparePointEstimate_y=nansum(comparePointEstimate_y,3);
comparePointEstimate_y(comparePointEstimate_y==0)=nan;

comparePointEstimate_z=cell2mat(permute({data.comparePointEstimate_z},[1 3 2]));
comparePointEstimate_z=nansum(comparePointEstimate_z,3);
comparePointEstimate_z(comparePointEstimate_z==0)=nan;

xyzRefAll=cell2mat(permute({data.xyzRefAll},[1 3 2]));
xyzRefAll=nansum(xyzRefAll,3);
xyzRefAll(xyzRefAll==0)=nan;

    xmean=nanmean(comparePointEstimate_x)';
    xstd=nanstd(comparePointEstimate_x)';
    ymean=nanmean(comparePointEstimate_y)';
    ystd=nanstd(comparePointEstimate_y)';
    zmean=nanmean(comparePointEstimate_z)';
    zstd=nanstd(comparePointEstimate_z)';
    compareAllX(:,:,i)=comparePointEstimate_x;
    compareAllY(:,:,i)=comparePointEstimate_y;
    compareAllZ(:,:,i)=comparePointEstimate_z;
    xyzRefAll_zscore=xyzRefAll-[xmean ymean zmean];
    xyzRefAll_zscore=xyzRefAll_zscore./[xstd ystd zstd];
    zDistance=sqrt(sum(xyzRefAll_zscore.^2,2));
        
    newX=nanmedian(comparePointEstimate_x);
    newY=nanmedian(comparePointEstimate_y);
    newZ=nanmedian(comparePointEstimate_z);
    
    newXAll(i,:)=newX;
    newYAll(i,:)=newY;
    newZAll(i,:)=newZ;
    
    zScoreAll(i,:)=zDistance;
    oldXAll(i,:)=xyzRefAll(:,1);
    oldYAll(i,:)=xyzRefAll(:,2);
    oldZAll(i,:)=xyzRefAll(:,3);
end
end


%[~,~,colorAll]=ndgrid(1:200,1:1522,1:93);
%% make volume plot
VolumeAll=nan(size(oldXAll));
for iTime=1:length(pointStats2);
    P1=pointStats2(iTime);
    trackIdx=P1.trackIdx;
    presentPoints=find(~isnan(trackIdx) & trackIdx<=size(oldXAll,1));
    VolumeAll(trackIdx(presentPoints),iTime)=P1.Volume(presentPoints);
    
    

    
end




%%

%%
window=7;
[searchX,searchY,searchZ]=ndgrid(-window:window,-window:window,-window:window);
searchX=searchX(:);
searchY=searchY(:);
searchZ=searchZ(:);

ctrlPoints=[searchX,searchY,searchZ];
newXAll2=nan(size(newXAll));
newYAll2=nan(size(newXAll));
newZAll2=nan(size(newXAll));
detAll=nan(size(newXAll));

for iTime=1:length(pointStats2)
    %%
    display(['Starting ' num2str(iTime)]);
    if ~isempty(pointStats2(iTime).stackIdx)
    imageFile=[imageFolder filesep 'image' num2str(pointStats2(iTime).stackIdx,'%4.5d') '.tif'];
    
    currentImageStack=stackLoad(imageFile,128); %have to put in image stack size
    imSize=size(currentImageStack);
    newVec=nan(size(compareAllX,3),3);
    detTemp=nan(size(compareAllX,3),1);
    for pointIdx=1:size(compareAllX,3)
        currentPoint=[oldXAll(pointIdx,iTime),oldYAll(pointIdx,iTime),...
            oldZAll(pointIdx,iTime)];
    allPoints=[compareAllX(:,iTime,pointIdx),compareAllY(:,iTime,pointIdx),...
        compareAllZ(:,iTime,pointIdx)];
    
    allPoints(any(isnan(allPoints),2) | any(allPoints==0,2),:)=[];
    allPoints(any(bsxfun(@ge,allPoints,imSize),2),:)=[];
    allPoints(any(allPoints<1,2),:)=[];
    
    if ~isempty(allPoints);
    allPointsR=round(allPoints);
    
   
    allPointsIdx=sub2ind(imSize,allPointsR(:,1),allPointsR(:,2),allPointsR(:,3));
    pointW=currentImageStack(allPointsIdx);
    floorLvl=quantile(pointW,.4);
    pointW=pointW-floorLvl;
    pointW(pointW<0)=0;
    if any(pointW)
    newMean=pointW'*allPoints/sum(pointW);
    newCov=bsxfun(@minus, allPoints,newMean);
    
    newCov=bsxfun(@times,newCov,pointW)'*newCov;
    newCov=newCov/sum(pointW);
    newCov=newCov*10;
    E=sum(((newCov)\ctrlPoints')'.*ctrlPoints,2);
    P=exp(-E);
    ctrlPoints2=round(bsxfun(@plus, ctrlPoints,newMean));
    inImage=~any(bsxfun(@gt, ctrlPoints2,imSize)| bsxfun(@lt, ctrlPoints2,[1 1 1]),2);
    ctrlPoints2Idx=sub2ind(imSize,ctrlPoints2(inImage,1),ctrlPoints2(inImage,2),ctrlPoints2(inImage,3));
    P=P(inImage);
    intensityWeight=currentImageStack(ctrlPoints2Idx);
    intensityWeight=intensityWeight-floorLvl;
    intensityWeight(intensityWeight<0)=0;
    totalWeight=P.*intensityWeight;
    totalWeight=totalWeight./sum(totalWeight);
    newMean2=newMean+[sum(searchX(inImage).*totalWeight),...
        sum(searchY(inImage).*totalWeight),...
        sum(searchZ(inImage).*totalWeight)];
    newVec(pointIdx,:)=newMean2;
detTemp(pointIdx)=det(newCov/10);
    else
        newVec(pointIdx,:)=nan;
        detTemp(pointIdx)=nan;
    end
    
    end
    end
    detAll(:,iTime)=detTemp;
    newXAll2(:,iTime)=newVec(:,1);
    newYAll2(:,iTime)=newVec(:,2);
    newZAll2(:,iTime)=newVec(:,3);
    end
    end

%% make comparison distance matrices
nanmat=isnan(newXAll2);
limit1=find(sum(~nanmat,2),1,'last');
limit2=find(sum(~nanmat,1),1,'last');

newXAll2=newXAll2(1:limit1,1:limit2);
newYAll2=newYAll2(1:limit1,1:limit2);
newZAll2=newZAll2(1:limit1,1:limit2);
detAll=detAll(1:limit1,1:limit2);
zScoreAll=zScoreAll(1:limit1,1:limit2);
for iTime=1:limit2
    P=[newXAll2(:,iTime), newYAll2(:,iTime), newZAll2(:,iTime)];
    Pold=[oldXAll(:,iTime), oldYAll(:,iTime), oldZAll(:,iTime)];
    dMat=squareform(pdist(P));
    dMatAll(:,:,iTime)=dMat;
    dMatOld=squareform(pdist(Pold));
    dMatAllOld(:,:,iTime)=dMat;
    
end

%% calculate zscores of distance matrices for new points, take average of highest 10
zMatAll=bsxfun(@minus, dMatAll,nanmean(dMatAll,3));
zMatSTD=nanstd(zMatAll,[],3);
zMatAll=bsxfun(@rdivide,zMatAll, zMatSTD);

zMatAll=sort(zMatAll,1,'descend');
zMatAll=zMatAll(1:10,:,:);
zMatAll=squeeze(nanmean(zMatAll,1));
imagesc(zMatAll)

zMatAllOld=bsxfun(@minus, dMatAllOld,nanmean(dMatAllOld,3));
zMatSTDOld=nanstd(zMatAllOld,[],3);
zMatAllOld=bsxfun(@rdivide,zMatAllOld, zMatSTDOld);

zMatAllOld=sort(zMatAllOld,1,'descend');
zMatAllOld=zMatAllOld(1:10,:,:);
zMatAllOld=squeeze(nanmean(zMatAllOld,1));
figure
imagesc(zMatAllOld)



%%
logDet=log(detAll);

check=(zScoreAll>3) & (logDet>6);
check(1:limit1,1:limit2)=check & zMatAll>6;
for i=1:length(pointStatsNew);
    iTime=i;%pointStatsNew(i).stackIdx;
    if any(check(:,iTime))
        replaceIdx=find(check(:,iTime));
        
        for iReplace=replaceIdx'
            lookupIdx=pointStatsNew(i).trackIdx==iReplace;
            if any(lookupIdx)
                pointStatsNew(i).straightPoints(lookupIdx,:)=...
[newXAll2(iReplace,iTime),newYAll2(iReplace,iTime),newZAll2(iReplace,iTime)];
                display([' Success? at ' num2str(iTime) ' ' num2str(iReplace)]);
            else
                display([' mismatch? at ' num2str(iTime) ' ' num2str(iReplace)]);
            end
            
        end
    end
        newIdx=find(isnan(zScoreAll(:,iTime)));
        newIdx=newIdx(~ismember(newIdx,pointStatsNew(i).trackIdx));
        for iNew=newIdx'
            lookupIdx=pointStatsNew(i).trackIdx==iNew;
            if ~any(lookupIdx)
          pointStatsNew(i).straightPoints=[  pointStatsNew(iTime).straightPoints;...
               [newXAll(iNew,iTime),newYAll(iNew,iTime),newZAll(iNew,iTime)]];
           pointStatsNew(i).trackIdx=[pointStatsNew(i).trackIdx;...
               iNew];
           display('New Point');
            elseif sum(lookupIdx)==1
        pointStatsNew(i).straightPoints(lookupIdx,:)=...
                    [newXAll2(iNew,iTime),newYAll2(iNew,iTime),newZAll2(iNew,iTime)];
           display('Replace Point');
            end   
        
        
        
        
    end
end


%% make average mindist matrix

dminAll=nan(size(newXAll2));

for iTime=1:size(newXAll2,2);
   currentData=[newXAll2(:,iTime) newYAll2(:,iTime) newZAll2(:,iTime)];
   D=squareform(pdist((currentData)));
   D(eye(size(D))>0)=nan;
   dmin=nanmean(D);
   
   dminAll(:,iTime)=dmin;
   
    
end

dmin=nanmean(dminAll,2);



%%

newXAll(newXAll==0)=nan;
newYAll(newXAll==0)=nan;
newZAll(newZAll==0)=nan;
% 
% for i=1:length(pointStatsNew);
%     iTime=i;% pointStatsNew(i).
%     pointStatsNew(i).straightPoints=[newXAll(:,i),newYAll(:,i),newZAll(:,i)];
% pointStatsNew(i).trackIdx=1:size(newXAll,1);
% %pointStatsNew(i).stackIdx=i;
% 
% end
% 

%% create pointStatsfor just one point

pointStatFolder=[dataFolder filesep 'PointStatFolder'];
mkdir(pointStatFolder);
for pointIdx=1:size(compareAllX,3);
P1stats=pointStats2;
display(['Starting ' num2str(pointIdx)]);
parfor iTime=1:length(P1stats)
    trackIdx=pointStats2(iTime).trackIdx;
    allPoints=[compareAllX(:,iTime,pointIdx),compareAllY(:,iTime,pointIdx),...
        compareAllZ(:,iTime,pointIdx)];
    idxList=nan(length(allPoints),1);
    if any(trackIdx==pointIdx);
        allPoints=[allPoints;  pointStats2(iTime).straightPoints(trackIdx==pointIdx,:)];
        idxList=[idxList; pointIdx];
    end
        trackIdx=pointStatsNew(iTime).trackIdx;

    allPoints=[allPoints;  pointStatsNew(iTime).straightPoints(trackIdx==pointIdx,:)];
    idxList=[idxList; 1];

    
    P1stats(iTime).straightPoints=allPoints;
      P1stats(iTime).trackIdx=idxList;
    P1stats(iTime).DMatrixi_x=[];
    P1stats(iTime).DMatrixi_y=[];
    P1stats(iTime).DMatrixi_z=[];
    
end


save([pointStatFolder filesep 'pointStatsTemp' num2str(pointIdx)],'P1stats');

end
%%
save([pointStatFolder filesep 'pointStatsNew'],'pointStatsNew');
