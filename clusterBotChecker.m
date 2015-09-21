function clusterBotChecker(filePath,startIdx,groupSize)
 %cluster bot checker uses a pointstats and a track index and stepsize as
 %inputs. The code goes through all frames and compares them to nSubSample
 %number of other randomly selected frames. In each frame, it asks where
 %the otherframes would guess where a certain neuron should be using TPS
 %interp.
 
 % changed to run under the 1 hr. 
 
 
if nargin==0
    filePath=uipickfiles;
    startIdx=1;
    filePath=filePath{1};
end
if nargin<3
    groupSize=0;
end
         outputName=fileparts(filePath);

load(filePath);
%%
%select number of points, 
nTime=length(pointStats2);
nSubSample=nTime;

startIdxReal=ceil(startIdx/groupSize);
runIdx=mod(startIdx,groupSize);
timeVector=1:nTime;
timeIdx=floor(startIdx/max(startIdx+1)*groupSize);
timeVector=timeVector(timeIdx==runIdx);
%How many neurons to check

% loop through all neurons in list
for iPointIdx=startIdxReal;
    display([' Starting ' num2str(iPointIdx)]);

comparePointEstimate_x=nan(nSubSample,nTime);
comparePointEstimate_y=nan(nSubSample,nTime);
comparePointEstimate_z=nan(nSubSample,nTime);
xyzRefAll=nan(nTime,3);

%loop through all time points
for iFrame=1:length(timeVector)
    %take random subsample of other time points to compare
   % subSample=randperm(nTime,nSubSample);
    subSample=1:nTime;
pointsRef=pointStats2(timeVector(iFrame));
refTrackIdx=pointsRef.trackIdx;
if any(refTrackIdx==iPointIdx)
xyzRef=pointsRef.straightPoints(refTrackIdx==iPointIdx,:);        
xyzRefAll(iFrame,:)=xyzRef;

end
%loop through comparison time points
for iCounter=1:nSubSample
    
    pointsI=pointStats2(subSample(iCounter));
    iTrackIdx=pointsI.trackIdx;
[~,overlapI,overlapRef]=intersect(iTrackIdx,refTrackIdx);

overlapI_excludeI=overlapI(iTrackIdx(overlapI)~=iPointIdx);
overlapRef_excludeI=overlapRef(refTrackIdx(overlapRef)~=iPointIdx);

    if any((overlapI)) && any(iTrackIdx==iPointIdx)
%get points excluding the point in question
        movingPoints_excludeI=pointsI.straightPoints(overlapI_excludeI,:);
        controlPoint_excludeI=pointsRef.straightPoints(overlapRef_excludeI,:);
xyzI=pointsI.straightPoints(iTrackIdx==iPointIdx,:);
%estimate where the moving frame thinks the reference frame's point should
%be

newEstimatePoint=tpswarp3points(movingPoints_excludeI,controlPoint_excludeI,xyzI);


comparePointEstimate_x(iCounter,iFrame)=newEstimatePoint(1);
comparePointEstimate_y(iCounter,iFrame)=newEstimatePoint(2);
comparePointEstimate_z(iCounter,iFrame)=newEstimatePoint(3);
    end

end

end
    outputNameFile=[outputName filesep 'botChecker' num2str(iPointIdx,'%3.5d')...
        'Run' num2str(runIdx,'%3.2d')];
    save(outputNameFile,'comparePointEstimate_x','comparePointEstimate_y',...
        'comparePointEstimate_z','xyzRefAll');
    
    
end