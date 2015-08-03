function clusterBotChecker(filePath,startIdx,stepSize)
 
if nargin==0
    filePath=uipickfiles;
    startIdx=1;
    filePath=filePath{1};
end
if nargin<3
    stepSize=0;
end
load(filePath);
%%
nSubSample=200;
nTime=length(pointStats2);

iIdxList=startIdx:(startIdx+stepSize);%:startIdx+stepSize-1;

for iPointIdx=iIdxList;
    display([' Starting ' num2str(iPointIdx)]);

comparePointEstimate_x=nan(nSubSample,nTime);
comparePointEstimate_y=nan(nSubSample,nTime);
comparePointEstimate_z=nan(nSubSample,nTime);
xyzRefAll=nan(nTime,3);
for iFrame=1:nTime
    subSample=randperm(nTime,nSubSample);
pointsRef=pointStats2(iFrame);
refTrackIdx=pointsRef.trackIdx;
if any(refTrackIdx==iPointIdx)
xyzRef=pointsRef.straightPoints(refTrackIdx==iPointIdx,:);        
xyzRefAll(iFrame,:)=xyzRef;

end

for iCounter=1:nSubSample
    
    pointsI=pointStats2(subSample(iCounter));
    iTrackIdx=pointsI.trackIdx;
[~,overlapI,overlapRef]=intersect(iTrackIdx,refTrackIdx);

overlapI_excludeI=overlapI(iTrackIdx(overlapI)~=iPointIdx);
overlapRef_excludeI=overlapRef(refTrackIdx(overlapRef)~=iPointIdx);

    if any((overlapI)) && any(iTrackIdx==iPointIdx)

     %       controlPoint=pointsRef.straightPoints(overlapRef,:);
      %      movingPointstemp=pointsI.straightPoints(overlapI,:);
    %    d2controlPoint=pdist2(controlPoint,movingPointstemp);
      %  d2controlPoint=d2controlPoint(tril(true(overlapN),-1))';

        movingPoints_excludeI=pointsI.straightPoints(overlapI_excludeI,:);
        controlPoint_excludeI=pointsRef.straightPoints(overlapRef_excludeI,:);
xyzI=pointsI.straightPoints(iTrackIdx==iPointIdx,:);
        
%tform = makeAffine3d(movingPoints_excludeI, controlPoint_excludeI);

% movingPoints_excludeI...
%     =transformPointsForward(tform,movingPoints_excludeI);
% 
% xyzI=transformPointsForward(tform,xyzI);

newEstimatePoint=tpswarp3points(movingPoints_excludeI,controlPoint_excludeI,xyzI);

comparePointEstimate_x(iCounter,iFrame)=newEstimatePoint(1);
comparePointEstimate_y(iCounter,iFrame)=newEstimatePoint(2);
comparePointEstimate_z(iCounter,iFrame)=newEstimatePoint(3);
    end

end

end
         outputName=fileparts(filePath);
    outputName=[outputName filesep 'botChecker' num2str(iPointIdx,'%3.5d')];
    save(outputName,'comparePointEstimate_x','comparePointEstimate_y',...
        'comparePointEstimate_z','xyzRefAll');
    
    
end