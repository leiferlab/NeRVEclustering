function clusterWormTracker(filePath,startIdx,nGroups,offset,doGroups)
%made specifically for 1hr queue, can only do ~ 250 comparisons per hour
 show=00;
 startTic=tic;
if nargin==0
    filePath=uipickfiles;
    startIdx=1;
    filePath=filePath{1};
end

if nargin<3
    nGroups=1;
end
if nargin<4
    offset=0;
end
if nargin<5
    doGroups=1;
end
if doGroups==1
 timeLimit=3600; 
else
    timeLimit=3600*doGroups*2;
end

startIdx=startIdx+offset;
load(filePath);
%%
 matchesPerSegment=150;
matchesPerSegment=matchesPerSegment*nGroups;
startIdx=(1:doGroups)+(startIdx-1)*doGroups;
iIdxList=floor((nGroups+startIdx-1)/nGroups);
itIdx=mod(startIdx,nGroups);
runIdxListAll=find(cellfun(@(x) ~isempty(x),{pointStats.stackIdx}));
presentN=length(runIdxListAll);
deltaRun=presentN/matchesPerSegment;


% presentIdx=cellfun(@(x) ~isempty(x),{pointStats.stackIdx},'uniform',0);
% presentIdx=find(cell2mat(presentIdx));
presentIdx=1:length(pointStats);
N=length(presentIdx);
param.dim=3;
param.good=2;
param.excessive=4;
param.quiet=1;
param.timeLimit=10;
param.difficult=1.5e4;
%%
for iCounter=1:length(iIdxList)%length(TrackData)
    %%
    iIdx=iIdxList(iCounter);
    runIdxList=runIdxListAll(ceil((deltaRun:deltaRun:presentN/nGroups)+itIdx(iCounter)*presentN/nGroups));
runIdxList=unique(runIdxList);

                outputName=fileparts(filePath);
                if isempty(outputName)
                    
    outputName=[outputName filesep 'trackMatrix' num2str(iIdx,'%3.5d') 'Run' num2str(itIdx(iCounter),'%3.2d')];
      display(outputName);

    i=presentIdx(iIdx);
    outRange=1:N;%max(1,i-windowSearch):min(length(TrackData),i+windowSearch);
    TrackMatrixi=zeros(size(pointStats(i).straightPoints,1),length(runIdxList));
    DMatrixi_x=TrackMatrixi;
    DMatrixi_y=TrackMatrixi;
    DMatrixi_z=TrackMatrixi;
    
    for runIdx=1:length(runIdxList)%outRange;
        %%
        j=runIdxList(runIdx);

        itic=tic;
        try
P1=pointStats(i);
P2=pointStats(j);
T1=[P1.straightPoints P1.Volume.^(1/3) P1.Rintensities];
T2=[P2.straightPoints P2.Volume.^(1/3) P2.Rintensities];

T1temp=pointStats(i).straightPoints(:,1:3);
T2temp=pointStats(j).straightPoints(:,1:3);

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = ...
    gmmreg_L2_multilevel_jn(T2...
    ,T1,1, [ 1,.5], ...
    [.000008, 0.0000008, 0.0008],[0 0],...
    [ 0.000001 0.0001 0.001 0.001],show);
trackInput=[T1temp T1temp  (1:length(T1temp))'  ones(size(T1temp(:,1))); ...
    Transformed_M(:,1:3) T2temp  (1:length(T2temp))' 2*ones(size(Transformed_M(:,1)))];
TrackOut=nan;
idx = kmeans(trackInput(:,1:3),3);
idx1=idx(1:length(T1temp));
idx2=idx(length(T1temp)+1:end);
%%

for iRegions=1:max(idx)

    [Transformed_M(idx2==iRegions,:), ~, multilevel_param] = ...
    gmmreg_L2_multilevel_jn(...
    Transformed_M(idx2==iRegions,:),T1(idx1==iRegions,:),  2, [ 3,.3,.3], ...
    [0.005,.0005, 0.002, 0.08],[0 0],...
    [ 0.000001 0.000001 0.000001 0.001],show);




end

%%

track1=[];
track2=[];

trackInput=[T1temp T1temp  (1:length(T1temp))'  ones(size(T1temp(:,1))); ...
    Transformed_M(:,1:3) T2temp  (1:length(T2temp))' 2*ones(size(Transformed_M(:,1)))];
TrackOut=nan;


for iRegions=1:max(idx)
    if any(idx==iRegions)
    trackInputi=trackInput(idx==iRegions,:);
 %   trackInputi(:,3)=trackInputi(:,3)*0;
    counter=18;
    TrackOut=nan;
    if length( unique(trackInputi(:,end)))>1
    while(all(isnan(TrackOut(:))))
        TrackOut=trackJN(trackInputi,counter,param);
        counter=counter-1;
    end
        
    TrackOut(:,1:3)=[];
    TrackStats=round(TrackOut(:,4:end));
   % TrackedIDs=TrackStats([1;diff(TrackStats(:,3))]==0,end);
    
  %  TrackStats=TrackStats(ismember(TrackStats(:,end),TrackedIDs),:);
    track1i=TrackStats(1:2:end,1);
    track2i=TrackStats(2:2:end,1);
    track1=[track1;track1i];
    track2=[track2;track2i];
    end
    end
end            
            
                    TrackMatrixi(track1,runIdx-outRange(1)+1)=track2;
                 matchedIdx=find( TrackMatrixi(:,runIdx-outRange(1)+1));
                matchedPairs=[matchedIdx, TrackMatrixi(matchedIdx,runIdx-outRange(1)+1)];
                
                
presentIJ=TrackMatrixi(:,runIdx-outRange(1)+1)>0;
points1=T1temp(track1,1:3);
points2=Transformed_M(track2,1:3);
pointsDiff=abs(points1-points2);

DMatrixi_x(presentIJ,runIdx-outRange(1)+1)=pointsDiff(:,1);
DMatrixi_y(presentIJ,runIdx-outRange(1)+1)=pointsDiff(:,2);
DMatrixi_z(presentIJ,runIdx-outRange(1)+1)=pointsDiff(:,3);

% 
% figure;
% scat3(T1temp);
% hold on
% scat3(Transformed_M);
% % track1=matchedPairs(:,1);
% % track2=matchedPairs(:,2);
% 
% plot3([T1temp(track1,1),Transformed_M(track2,1)]',...
%     [T1temp(track1,2),Transformed_M(track2,2)]',...
%     [T1temp(track1,3),Transformed_M(track2,3)]','linewidth',4)
% axis equal
% 
% 
% figure
% scat3(T1temp);
% hold on
% scat3(T2temp);
% % track1=matchedPairs(:,1);
% % track2=matchedPairs(:,2);
% 
% plot3([T1temp(track1,1),T2temp(track2,1)]',...
%     [T1temp(track1,2),T2temp(track2,2)]',...
%     [T1temp(track1,3),T2temp(track2,3)]','linewidth',4)
% axis equal

      display(['Finished match' num2str(j) ' in ' num2str((toc(itic))) 's']);

      %for the hr queue on cluster, , if run time is greater than an hour,
      %time to save at each iteration. 

if toc(startTic)>timeLimit;
        save(outputName,'TrackMatrixi');    
end
        catch ME
            ME
        end
    end
    if isempty(TrackMatrixi)
        TrackMatrixi=[];
    end
    
    save(outputName,'TrackMatrixi');
end

