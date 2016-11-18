function clusterBotCheckCompiler(dataFolder)
% clusterBotCheckCompiler takes the results of clusterBotChecker and and
% and uses the points to correct / fill in gaps in the pointStats file
% produced but wormTrackCompiler. The input is that dataFolder that has the
% PointStast2 file and the BotCheck folder with all of the BotCheck files.

%% load inputs
if nargin==0
    %if no inputs given, manually select them.
    display('Select submission Folder')
    mostRecent=getappdata(0,'mostRecent');
    display('Select Botchecker Folder');
    submissionFolder=uipickfiles('filterSpec',mostRecent);
    setappdata(0,'mostRecent',submissionFolder{1});
    submissionFolder=submissionFolder{1};
    display('Select pointstats')
    pointStatsFile=uipickfiles('filterspec',fileparts(submissionFolder));
    pointStatsFile=pointStatsFile{1};
    imageFolder=uipickfiles('filterSpec',fileparts(submissionFolder));
    imageFolder=imageFolder{1};
    dataFolder=fileparts(imageFolder);
else
    submissionFolder=[dataFolder filesep 'botCheckFolder'];
    pointStatsFile=[dataFolder filesep 'PointsStats2'];
    
    psFolder=dir([dataFolder filesep 'CLstraight*']);
    imageFolder=[dataFolder filesep psFolder.name];
    display(imageFolder);
end

%get list of all botcheck files
fileListAll=dir([submissionFolder filesep 'bot*.mat']);

%load in pointStats file
load(pointStatsFile);
pointStatsNew=pointStats2;




%%
%get subsample of points for position guessing
nsubCompare=500;
%get number of neurons, as the string at the end of the last botcheck file
nPoint=str2double(fileListAll(end).name(11:15));
%load one of the datasets for initialization purposes
data=load([submissionFolder filesep fileListAll(10).name]);

%get number of
[nCompare,n_times]=size(data.comparePointEstimate_x);

%initialize
oldXAll=nan(nPoint,n_times);
oldYAll=oldXAll;
oldZAll=oldXAll;
zScoreAll=oldXAll;
%initialize 3D matrix which will contain all the guesses for each point.
%the matrix is organized as comparison frame, current frame being checked,
%current point being checked.
compareAllX=zeros(nsubCompare,n_times,nPoint);
compareAllY=compareAllX;
compareAllZ=compareAllX;
weightAll=compareAllX;

subCompare=round(1:nCompare/nsubCompare:nCompare);
subCompare=unique(subCompare);
%% loop over all points and load all of the bot checker data
for i=1:nPoint
    %%
    
    display(['Loading file ' num2str(i)]);
    fileName=['botChecker' num2str(i,'%3.5d') '*'];
    %get file all files that check neuron i, accounting for the possibility
    %of multiple runs
    fileList=dir([submissionFolder filesep fileName]);
    fileList={fileList.name};
    % get number of runs for each nuron
    nChecks=length(fileList);
    %initialize entire matrix to house all guesses
    comparePointEstimate_x=nan(nCompare,nCompare,nChecks);
    comparePointEstimate_y=comparePointEstimate_x;
    comparePointEstimate_z=comparePointEstimate_x;
    comparePointsW=comparePointEstimate_x;
    xyzRefAll=nan(nCompare,3,nChecks);
    
    %load all of that data
    clear data
    if ~isempty(fileList)
        %load each run for a neuron
        for j=1:nChecks
            display(['loading ' fileName ' frame ' num2str(j) ]);
            data=load([submissionFolder filesep fileList{j}]);
            %load x, y, z and conf for each run, compiling them together in
            %the 3rd dimension. There shouldnt be overlap but this was done
            %just in case, I can probably make this better now.
            comparePointEstimate_x(:,:,j)=data.comparePointEstimate_x;
            comparePointEstimate_y(:,:,j)=data.comparePointEstimate_y;
            comparePointEstimate_z(:,:,j)=data.comparePointEstimate_z;
            comparePointsW(:,:,j)=data.comparePointConf;
            xyzRefAll(:,:,j)=data.xyzRefAll;
        end
        %take subset for each of these the comparison of each point to
        %every other point,
        comparePointEstimate_x=comparePointEstimate_x(subCompare,:,:);
        comparePointEstimate_y=comparePointEstimate_y(subCompare,:,:);
        comparePointEstimate_z=comparePointEstimate_z(subCompare,:,:);
        comparePointsW=comparePointsW(subCompare,:,:);
        
        %% set nans for each of these to zero and project to get single
        %matrix for
        nanmap=isnan(comparePointEstimate_x);
        comparePointEstimate_x(nanmap)=0;
        comparePointEstimate_x=sum(comparePointEstimate_x,3);
        
        comparePointEstimate_y(nanmap)=0;
        comparePointEstimate_y=sum(comparePointEstimate_y,3);
        
        comparePointEstimate_z(nanmap)=0;
        comparePointEstimate_z=sum(comparePointEstimate_z,3);
        
        comparePointsW(nanmap)=0;
        comparePointsW=sum(comparePointsW,3);
        
        xyzRefAll=nansum(xyzRefAll,3);
        xyzRefAll(xyzRefAll==0)=nan;
        
        comparePointEstimate_z(comparePointEstimate_z==0)=nan;
        comparePointEstimate_y(comparePointEstimate_y==0)=nan;
        comparePointEstimate_x(comparePointEstimate_x==0)=nan;
        
        %% find weighted average of coordinates
        
        % get weight normalization
        sumW=nansum(comparePointsW);
        %get weighted mean
        xmean=sum(comparePointEstimate_x.*comparePointsW)./sumW;
        ymean=sum(comparePointEstimate_y.*comparePointsW)./sumW;
        zmean=sum(comparePointEstimate_z.*comparePointsW)./sumW;
        %get weighted mean^2
        x2mean=sum(comparePointEstimate_x.^2.*comparePointsW)./sumW;
        y2mean=sum(comparePointEstimate_y.^2.*comparePointsW)./sumW;
        z2mean=sum(comparePointEstimate_z.^2.*comparePointsW)./sumW;
        %get variance
        xstd=sqrt(-xmean.^2+x2mean);
        ystd=sqrt(-ymean.^2+y2mean);
        zstd=sqrt(-zmean.^2+z2mean);
        
        %find zscore for each point relative to the point cloud of guesses
        xyzRefAll_zscore=xyzRefAll-[xmean' ymean' zmean'];
        xyzRefAll_zscore=xyzRefAll_zscore./[xstd' ystd' zstd'];
        %find mahalanobis distance of each found to the center of the point
        %cloud
        zDistance=sqrt(sum(xyzRefAll_zscore.^2,2));

        %% compile all results

        % compile mahalanobis distance of each point from the point cloud
        zScoreAll(i,:)=zDistance;
        
        %compile all old points
        oldXAll(i,:)=xyzRefAll(:,1);
        oldYAll(i,:)=xyzRefAll(:,2);
        oldZAll(i,:)=xyzRefAll(:,3);
        
        %compile all weights
        weightAll(:,:,i)=comparePointsW;
        
        % compareAllX(a,b,c) estiamtes the X position of point c in frame
        % a based on frame b
        compareAllX(:,:,i)=comparePointEstimate_x;
        compareAllY(:,:,i)=comparePointEstimate_y;
        compareAllZ(:,:,i)=comparePointEstimate_z;
    end
end

%% make grid of control points
box=7;
[searchX,searchY,searchZ]=ndgrid(-box:box,-box:box,-box:box);
searchX=searchX(:);
searchY=searchY(:);
searchZ=searchZ(:);
%%
pointStatsNew=pointStats2;
ctrlPoints=[searchX,searchY,searchZ];
newXAll2=nan(size(oldXAll));
newYAll2=nan(size(oldXAll));
newZAll2=nan(size(oldXAll));
detAll=nan(size(oldXAll));
%%
empty_frames=cellfun(@(x) ~isempty(x), {pointStats2.straightPoints});
empty_frames=find(empty_frames);
for iTime=empty_frames
    try
    %% loop over time points
    display(['Starting ' num2str(iTime)]);
    %    save([dataFolder filesep 'errorCatch'],'iTime');
    %    debug, trying to find
    %    out why catostrophic failure occurs
    stackIdx_str=num2str(pointStats2(iTime).stackIdx,'%4.5d');
    imageFile=[imageFolder filesep 'image' stackIdx_str '.tif'];
    psFile=[imageFolder filesep 'pointStats' stackIdx_str '.mat'];
    % load image and lookup tables
    %TODO: have to put in image stack size, will find way to not hardcode
    currentImageStack=stackLoad(imageFile,128);
    currentImageStack=normalizeRange(currentImageStack);
    imSize=size(currentImageStack);
    
    %load pointstats file corresponding to this stack
    subStackData=load(psFile);
    subStackData=subStackData.pointStats;
    %get the binary mask
    stackMask=subStackData.baseImg;
    currentImageStack=currentImageStack+stackMask*2;
    
    %get current pointstats (load from ps file is not updated with
    %track idx)
    currentPS=pointStats2(iTime);
    %get unannotated points for possible reassignment, getting the point
    %coordinates and their index
    track_idx=currentPS.trackIdx;
    unAnnotatedPoints=currentPS.straightPoints(isnan(track_idx),[1 2 3]);
    unAnnotatedIdx=find(isnan(track_idx));
    newVec=nan(size(compareAllX,3),3);
    detTemp=nan(size(compareAllX,3),1);
    %loop over points inside a given time
    for pointIdx=1:size(compareAllX,3)
        %all the guesses for a given point, at a given time
        allPoints=[compareAllX(:,iTime,pointIdx),...
            compareAllY(:,iTime,pointIdx),...
            compareAllZ(:,iTime,pointIdx)];
        %all the weights for the guesses
        weights=weightAll(:,iTime,pointIdx,:);
        
        %remove points and weights with nans, and points outside the range
        outside=bsxfun(@ge,allPoints,imSize)| allPoints==0 |allPoints<1;
        nan_points=isnan(allPoints) ;
        removes=outside | nan_points;
        removes=any(removes,2) | isnan(weights);
        %delete those points
        allPoints(removes,:)=[];
        weights(removes)=[];
        
        %find image intensity at all guessed points
        allPointsR=round(allPoints);
        allPointsIdx=sub2ind(...
            imSize,allPointsR(:,1),allPointsR(:,2),allPointsR(:,3));
        % weights is image intensity times fit weight
        pointW=(currentImageStack(allPointsIdx).^2).*weights;
        %ignore lower 40% of weights
        floorLvl=quantile(pointW,.4);
        pointW=pointW-floorLvl;
        pointW(pointW<0)=0;
        
        if any(pointW)
            %find weighted centroid
            newMean=pointW'*allPoints/sum(pointW);
            
            %get covariance of point cloud
            newCov=bsxfun(@minus, allPoints,newMean);
            newCov=bsxfun(@times,newCov,pointW)'*newCov;
            newCov=newCov/sum(pointW);
            newCov=newCov*10;
            
            %get mahalobonis distance of each unannotated point, if there
            %is a close one, assign it to that cluster
            mahD=bsxfun(@minus,unAnnotatedPoints,newMean);
            mahD=sqrt(sum(mahD/newCov.*mahD,2));
            
            %using the covariance of the guess cloud, calculate the closes 
            %point in terms of mahabalonis distance
            [mahD, closestPoint]=min(mahD);
            %if a point is closer than 1.5 cov away from the mean of the 
            % guesses but previously unassigned, add the assignment.
            closestPointIdx=unAnnotatedIdx(closestPoint);
            
            %if any of the points are closer than 1.5 std and are not
            %annotated somewhere else in the volume, set that point to be a
            %member of that cluster.
            if any(mahD<1.5) && ~any(track_idx==pointIdx)
                currentPS.trackIdx(closestPointIdx)=pointIdx;
                currentPS.trackWeights(closestPointIdx)=.5;
                pointStatsNew(iTime)=currentPS;
                newVec(pointIdx,:)=nan;
                detTemp(pointIdx)=nan;
            else
                %if none of the points are close, we will have to put down
                %a new point as the weighted centroid of the control points
                %weighted by both image intensity and the ellipsoidal
                %gaussian from the point cloud
                
                % Energy of each control point as exp(-r^2)
                P=exp(-sum((newCov\ctrlPoints')'.*ctrlPoints,2));
                
                %translate control grid to center around new mean
                ctrlPoints2=round(bsxfun(@plus, ctrlPoints,newMean));
                %get all control points in the image
                outside=bsxfun(@gt, ctrlPoints2,imSize);
                outside=outside| bsxfun(@lt, ctrlPoints2,[1 1 1]);
                inImage=~any(outside,2);
                ctrlPoints2=ctrlPoints2(inImage,:);
                ctrlPoints2Idx=sub2ind(imSize,...
                    ctrlPoints2(:,1),ctrlPoints2(:,2),ctrlPoints2(:,3));
                P=P(inImage);
                
                %get image intensities at all of control points
                intensityWeight=currentImageStack(ctrlPoints2Idx);
                %subtract a floor
                floorLvl=quantile(intensityWeight(:),.4);
                intensityWeight=intensityWeight-floorLvl;
                intensityWeight(intensityWeight<0)=0;
                
                %combine weights from point cloud and image intensity
                totalWeight=P.*intensityWeight;
                totalWeight=totalWeight./sum(totalWeight);
                
                %calculate new weighted mean by adding the change from new
                %weighting to the new centroid of the point cloud. 
                newMean2=newMean+[sum(searchX(inImage).*totalWeight),...
                    sum(searchY(inImage).*totalWeight),...
                    sum(searchZ(inImage).*totalWeight)];
                
                %compile results for point coordinates and the determinant
                %of the covariance matrix. 
                newVec(pointIdx,:)=newMean2;
                detTemp(pointIdx)=det(newCov/10);
            end
        else
            newVec(pointIdx,:)=nan;
            detTemp(pointIdx)=nan;
        end
        
    end
    %compile results for point coordinates and the determinant
    %of the covariance matrix for all points at all times
    detAll(:,iTime)=detTemp;
    newXAll2(:,iTime)=newVec(:,1);
    newYAll2(:,iTime)=newVec(:,2);
    newZAll2(:,iTime)=newVec(:,3);
    catch me
        me
    end
        
end

%% make comparison distance matrices, for all times, calculate the d matrix

for iTime=1:length(pointStatsNew) %limit2
    %get all pairwise distances between points in a frame
    P=[newXAll2(:,iTime), newYAll2(:,iTime), newZAll2(:,iTime)];
    dMat=squareform(pdist(P));
    dMatAll(:,:,iTime)=dMat;
    
end

%% calculate zscores of distance matrices for new points, 
% look over time of all the pairwise distances, and for each one, get a
% zScore
zMatAll=bsxfun(@minus,dMatAll,nanmean(dMatAll,3));
zMatSTD=nanstd(zMatAll,[],3);
zMatAll=bsxfun(@rdivide,zMatAll, zMatSTD);

%set the nans to -Inf
zMatAll(isnan(zMatAll))=-Inf;

%sort the pairwise distance zscores for each point in each frame
zMatAll=sort(zMatAll,1,'descend');

%% try to quantify if a point is bad
%for each time point, and for each neuron, look at the top 15 distances
%between that neuron and another neuron (by zScore). If the average of
%those zscores is large, something is suspcious with that point and it
%should be checked. 

%if a point has gone far off where it should
%be, then many of the pairwise distances will be off, if a different point is
%mislabelled, only one of the distances will be off, so the trim mean of
%the top 15 should ignore that. 

zMatAll=zMatAll(1:15,:,:);
%take the average of those top 15 points, for each point, get the top 15 z
%scores of the distances from that point to 
zMatAll=squeeze(trimmean(zMatAll,10,1));


%% check point cloud
%if the point cloud (indicated by the log determanant of the point cloud) is
%too large or if the point is far from the point cloud, check the point
logDet=log(abs(detAll));


check_pc=(zScoreAll>3) & (logDet>6);

%combine point cloud check with zmat check
check=(check_pc & zMatAll>3)| zMatAll>6;
%% replace terms if check criteria met and fill inblanks
if ~exist('pointStats3','var')
    pointStats3=pointStatsNew;
end
nanmat=isnan(newXAll2);
limit1=find(sum(~nanmat,2),1,'last');

% add in new points
for iTime=1:length(pointStatsNew);
    trackIdx=pointStatsNew(iTime).trackIdx;
    %add if there are any checks in that time
    if any(check(:,iTime))
        %find indices that need to be replaced
        replaceIdx=find(check(:,iTime));
        for iReplace=replaceIdx'
            %get the in volume neuron to be replaced
            lookupIdx=trackIdx==iReplace;
            if any(lookupIdx)
                %get the new point and new weight, then plug them in
                new_point=[...
                    newXAll2(iReplace,iTime),...
                    newYAll2(iReplace,iTime),...
                    newZAll2(iReplace,iTime)...
                    ];
                
                new_weight=-detAll(iReplace,iTime);
                pointStatsNew(iTime).straightPoints(lookupIdx,:)=new_point;
                pointStatsNew(iTime).trackWeights(lookupIdx)=new_weight;
                display(['success at ' num2str(iTime) '' num2str(iReplace)]);
            else
                display(['error at ' num2str(iTime) '' num2str(iReplace)]);
            end
        end
    end
    %new indexes are ones that are nans and not currently present (they're
    %missing and need to be re added)
    newIdx=find(isnan(zScoreAll(1:limit1,iTime)));
    newIdx=newIdx(~ismember(newIdx,trackIdx));
    for iNew=newIdx'
        lookupIdx=trackIdx==iNew;
        % if the current trackIdx  is missing, a new point needs to be
        % added at the end of the straight points and trackIdx
        if ~any(lookupIdx)
            display('New Point');
            new_point=[...
                newXAll2(iNew,iTime),...
                newYAll2(iNew,iTime),...
                newZAll2(iNew,iTime)...
                ];
            %add point to its pointStats
            pointStatsNew(iTime).straightPoints=[...
                pointStatsNew(iTime).straightPoints;...
                new_point...
                ];
            
            pointStatsNew(iTime).trackIdx=[...
                pointStatsNew(iTime).trackIdx;...
                iNew];
            
            pointStatsNew(iTime).trackWeights(lookupIdx)=-detAll(iNew,iTime);
        elseif sum(lookupIdx)==1
            display('Replace Point');
            %if the current trackIdx is not missing, just move the point
            new_point=[...
                newXAll2(iNew,iTime),...
                newYAll2(iNew,iTime),...
                newZAll2(iNew,iTime)...
                ];
            
            pointStatsNew(iTime).straightPoints(lookupIdx,:)=new_point;
            pointStatsNew(iTime).trackWeights(lookupIdx)=-detAll(iNew,iTime);
        end
    end
end

%% save results
save([dataFolder filesep 'pointStatsNew'],'pointStatsNew');

