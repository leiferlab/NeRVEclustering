display('Select submission Folder')
mostRecent=getappdata(0,'mostRecent');
display('Select Botchecker Folder');
submissionFolder=uipickfiles('filterSpec',mostRecent);
setappdata(0,'mostRecent',submissionFolder{1});
submissionFolder=submissionFolder{1};
fileListAll=dir([submissionFolder filesep 'bot*.mat']);
display('Select pointstats')
pointStatsFile=uipickfiles('filterspec',fileparts(submissionFolder));
%
load(pointStatsFile{1});
pointStatsNew=pointStats2;

% select image folder
display('Select straightened data folder');

imageFolder=uipickfiles('filterSpec',fileparts(submissionFolder));
imageFolder=imageFolder{1};
dataFolder=fileparts(imageFolder);

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
      size(data.comparePointEstimate_x,2),iPoint);
  compareAllY=compareAllX;
  compareAllZ=compareAllX;
    weightAll=compareAllX;

  nCompare=size(compareAllX,1);
  %%
for i=1:iPoint
    %%
   % i=pointStats2(iPS).stackIdx;
    fileName=['botChecker' num2str(i,'%3.5d') '*'];
    fileList=dir([submissionFolder filesep fileName]);
    fileList={fileList.name};
    nChecks=length(fileList);
        comparePointEstimate_x=nan(nCompare,nCompare,nChecks);
    comparePointEstimate_y=comparePointEstimate_x;
    comparePointEstimate_z=comparePointEstimate_x;
    comparePointsW=comparePointEstimate_x;
    xyzRefAll=nan(nCompare,3,nChecks);
 clear data
 if ~isempty(fileList)
for j=1:nChecks
    display(['loading ' fileName ' frame ' num2str(j) ]);
    data=load([submissionFolder filesep fileList{j}]);
    comparePointEstimate_x(:,:,j)=data.comparePointEstimate_x;
    comparePointEstimate_y(:,:,j)=data.comparePointEstimate_y;
    comparePointEstimate_z(:,:,j)=data.comparePointEstimate_z;
    comparePointsW(:,:,j)=data.comparePointConf;
    xyzRefAll(:,:,j)=data.xyzRefAll;


end
nanmap=isnan(comparePointEstimate_x);
comparePointEstimate_x(nanmap)=0;
comparePointEstimate_x=sum(comparePointEstimate_x,3);

comparePointEstimate_y(nanmap)=0;
comparePointEstimate_y=sum(comparePointEstimate_y,3);

comparePointEstimate_z(nanmap)=0;
comparePointEstimate_z=sum(comparePointEstimate_z,3);

comparePointsW(nanmap)=0;
comparePointsW=sum(comparePointsW,3);

sumW=nansum(comparePointsW);
xyzRefAll=nansum(xyzRefAll,3);
xyzRefAll(xyzRefAll==0)=nan;

    xmean=sum(comparePointEstimate_x.*comparePointsW)./sumW;
    x2mean=sum(comparePointEstimate_x.^2.*comparePointsW)./sumW;
    xstd=sqrt(-xmean.^2+x2mean);
    ymean=sum(comparePointEstimate_y.*comparePointsW)./sumW;
    y2mean=sum(comparePointEstimate_y.^2.*comparePointsW)./sumW;
    ystd=sqrt(-ymean.^2+y2mean);
    zmean=sum(comparePointEstimate_z.*comparePointsW)./sumW;
    z2mean=sum(comparePointEstimate_z.^2.*comparePointsW)./sumW;
    zstd=sqrt(-zmean.^2+z2mean);
comparePointEstimate_z(comparePointEstimate_z==0)=nan;
comparePointEstimate_y(comparePointEstimate_y==0)=nan;
comparePointEstimate_x(comparePointEstimate_x==0)=nan;



compareAllX(:,:,i)=comparePointEstimate_x;
    compareAllY(:,:,i)=comparePointEstimate_y;
    compareAllZ(:,:,i)=comparePointEstimate_z;
    weightAll(:,:,i)=(comparePointsW);
    xyzRefAll_zscore=xyzRefAll-[xmean' ymean' zmean'];
    xyzRefAll_zscore=xyzRefAll_zscore./[xstd' ystd' zstd'];
    zDistance=sqrt(sum(xyzRefAll_zscore.^2,2));
        
    newX=nanmedian(comparePointEstimate_x);
    newY=nanmedian(comparePointEstimate_y);pointStats2
    
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
pointStatsNew=pointStats2;
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
    %% loop over time points
    display(['Starting ' num2str(iTime)]);
    save([dataFolder filesep 'errorCatch'],'iTime');
    if ~isempty(pointStats2(iTime).stackIdx)
    imageFile=[imageFolder filesep 'image' num2str(pointStats2(iTime).stackIdx,'%4.5d') '.tif'];
    psFile=[imageFolder filesep 'pointStats' num2str(pointStats2(iTime).stackIdx,'%4.5d') '.mat'];
    % load image and lookup tables
    currentImageStack=stackLoad(imageFile,128); %have to put in image stack size
    currentImageStack=normalizeRange(currentImageStack);
    
    subStackData=load(psFile);
    subStackData=subStackData.pointStats;
    stackMask=subStackData.baseImg;
    currentImageStack=currentImageStack+stackMask*2;
    currentPS=pointStats2(iTime);
    unAnnotatedPoints=currentPS.straightPoints(isnan(currentPS.trackIdx),[1 2 3]);
    unAnnotatedIdx=find(isnan(currentPS.trackIdx));
    imSize=size(currentImageStack);
    newVec=nan(size(compareAllX,3),3);
    detTemp=nan(size(compareAllX,3),1);
    %loop over points inside a given time
    for pointIdx=1:size(compareAllX,3)
       
        currentPoint=[oldXAll(pointIdx,iTime),oldYAll(pointIdx,iTime),...
            oldZAll(pointIdx,iTime)];
    allPoints=[compareAllX(:,iTime,pointIdx),compareAllY(:,iTime,pointIdx),...
        compareAllZ(:,iTime,pointIdx)]; %all the guesses for a given point
    weights=weightAll(:,iTime,pointIdx,:); %all the weights for the guesses
    %remove nans, and points outside the range
    removes=isnan(allPoints) | allPoints==0 |allPoints<1 |bsxfun(@ge,allPoints,imSize);
    removes=any(removes,2) | isnan(weights);
    allPoints(removes,:)=[];
    weights(removes)=[];
    
    if ~isempty(allPoints) %if there are any guesses
        %find image intensity at all guessed points
    allPointsR=round(allPoints);
    allPointsIdx=sub2ind(imSize,allPointsR(:,1),allPointsR(:,2),allPointsR(:,3));
   % weights is image intensity times fit weight
    pointW=(currentImageStack(allPointsIdx).^2).*weights;
    %ignore lower 40% of weights
        floorLvl=quantile(pointW,.4);
    pointW=pointW-floorLvl;
    pointW(pointW<0)=0;
    
    if any(pointW)
    newMean=pointW'*allPoints/sum(pointW);
    newCov=bsxfun(@minus, allPoints,newMean);
    
    newCov=bsxfun(@times,newCov,pointW)'*newCov;
    newCov=newCov/sum(pointW);
    newCov=newCov*10;
    d_currentPoints=bsxfun(@minus,unAnnotatedPoints,newMean);
    currentPoints_mahD=sqrt(sum(d_currentPoints/newCov .*d_currentPoints,2));
    %using the covariance of the guess cloud, calculate the closes point in
    %terms of mahabalonis distance
    [mahD, closestPoint]=min(currentPoints_mahD);
    closestPointIdx=unAnnotatedIdx(closestPoint);
    %if a point is closer than 1.5 cov away from the mean of the guesses
    %but previously unassigned, add the assignment.
    
    if any(mahD<1.5) && ~any(currentPS.trackIdx==pointIdx) 
currentPS.trackIdx(closestPointIdx)=pointIdx;
currentPS.trackIdx(currentPS.trackIdx==pointIdx)=NaN;


currentPS.trackWeights(closestPointIdx)=currentPS.trackWeights(pointIdx);
currentPS.trackWeights(pointIdx)=0;

pointStatsNew(iTime)=currentPS;

        newVec(pointIdx,:)=nan;
        detTemp(pointIdx)=nan;
    else
        
    E=sum(((newCov)\ctrlPoints')'.*ctrlPoints,2);
    P=exp(-E);
    ctrlPoints2=round(bsxfun(@plus, ctrlPoints,newMean));
    inImage=~any(bsxfun(@gt, ctrlPoints2,imSize)| bsxfun(@lt, ctrlPoints2,[1 1 1]),2);
    ctrlPoints2Idx=sub2ind(imSize,ctrlPoints2(inImage,1),ctrlPoints2(inImage,2),ctrlPoints2(inImage,3));
    P=P(inImage);
    intensityWeight=currentImageStack(ctrlPoints2Idx);
        floorLvl=quantile(intensityWeight(:),.4);

    intensityWeight=intensityWeight-floorLvl;
    intensityWeight(intensityWeight<0)=0;
    totalWeight=P.*intensityWeight;
    totalWeight=totalWeight./sum(totalWeight);
    newMean2=newMean+[sum(searchX(inImage).*totalWeight),...
        sum(searchY(inImage).*totalWeight),...
        sum(searchZ(inImage).*totalWeight)];
    
    newVec(pointIdx,:)=newMean2;
detTemp(pointIdx)=det(newCov/10);
    end
    else
        newVec(pointIdx,:)=nan;
        detTemp(pointIdx)=nan;
    end
    
    end
    end
    
     trackIdx=pointStatsNew(iTime).trackIdx;
    trackIdx=trackIdx(~isnan(trackIdx));
    if length(trackIdx)>length(unique(trackIdx))
        keyboard
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
% 
% newXAll2=newXAll2(1:limit1,1:limit2);
% newYAll2=newYAll2(1:limit1,1:limit2);
% newZAll2=newZAll2(1:limit1,1:limit2);
% detAll=detAll(1:limit1,1:limit2);
% zScoreAll=zScoreAll(1:limit1,1:limit2);
for iTime=1:length(pointStatsNew) %limit2
    P=[newXAll2(:,iTime), newYAll2(:,iTime), newZAll2(:,iTime)];
    Pold=[oldXAll(:,iTime), oldYAll(:,iTime), oldZAll(:,iTime)];
    dMat=squareform(pdist(P));
    dMatAll(:,:,iTime)=dMat;
    dMatOld=squareform(pdist(Pold));
    dMatAllOld(:,:,iTime)=dMat;
    
end

%% calculate zscores of distance matrices for new points, take average of highest 10
%check out the distnace matrices,convert to zscores
zMatAll=bsxfun(@minus, dMatAll,nanmean(dMatAll,3));
zMatSTD=nanstd(zMatAll,[],3);
zMatAll=bsxfun(@rdivide,zMatAll, zMatSTD);
zMatAll(isnan(zMatAll))=-Inf;
%sort the pairwise distance zscores for each point in each frame
zMatAll=sort(zMatAll,1,'descend');
%take only the top 15 zscores, if a point has gone far off where it should
%be, the many of the pairwise distances will be off, if another point is
%mislabelled, only one of the distsances will be off .
zMatAll=zMatAll(1:15,:,:);
zMatAll=squeeze(trimmean(zMatAll,10,1));



%%
logDet=log(abs(detAll));


check=(zScoreAll>3) & (logDet>6);
check=(check & zMatAll>3) | zMatAll>6;
%check(1:limit1,1:limit2)=(check & zMatAll>3) | zMatAll>6;

%% replace terms if check criteria met and fill inblanks
if ~exist('pointStats3','var')
pointStats3=pointStatsNew;
end

for i=1:length(pointStatsNew);
    iTime=i;%pointStatsNew(i).stackIdx;
    if any(check(:,iTime))
        replaceIdx=find(check(:,iTime));
        for iReplace=replaceIdx'
            lookupIdx=pointStatsNew(i).trackIdx==iReplace;
            if any(lookupIdx)
                pointStatsNew(i).straightPoints(lookupIdx,:)=...
[newXAll2(iReplace,iTime),newYAll2(iReplace,iTime),newZAll2(iReplace,iTime)];
                   pointStatsNew(i).trackWeights(lookupIdx)=-detAll(iReplace,iTime);
                display([' Success? at ' num2str(iTime) ' ' num2str(iReplace)]);
            else
                display([' mismatch? at ' num2str(iTime) ' ' num2str(iReplace)]);
            end
            
        end
    end
    % if point needs to be added beyond the number of points currently
    % present in the frame
        newIdx=find(isnan(zScoreAll(1:limit1,iTime)));
        newIdx=newIdx(~ismember(newIdx,pointStatsNew(i).trackIdx));
        for iNew=newIdx'
            lookupIdx=pointStatsNew(i).trackIdx==iNew;
            if ~any(lookupIdx)
          pointStatsNew(i).straightPoints=[  pointStatsNew(iTime).straightPoints;...
               [newXAll(iNew,iTime),newYAll(iNew,iTime),newZAll(iNew,iTime)]];
           pointStatsNew(i).trackIdx=[pointStatsNew(i).trackIdx;...
               iNew];
                   pointStatsNew(i).trackWeights(lookupIdx)=-detAll(iReplace,iTime);

           display('New Point');
            elseif sum(lookupIdx)==1
        pointStatsNew(i).straightPoints(lookupIdx,:)=...
                    [newXAll2(iNew,iTime),newYAll2(iNew,iTime),newZAll2(iNew,iTime)];
                   pointStatsNew(i).trackWeights(lookupIdx)=-detAll(iReplace,iTime);

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
save([dataFolder filesep 'pointStatsNew'],'pointStatsNew');


