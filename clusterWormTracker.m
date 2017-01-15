function clusterWormTracker(filePath,startIdx,nRef,start_offset)
% clusterWormTracker compares a set of pointsets from WormStraighten code
% and uses non rigid pointset registration to create match matrices. The
% code has been now modified to just run with on one sample with nRef for a
% given pointStats. You can no longer split up a single sample.

%%%% Inputs
% filePath : link to the complete PointStats file with data from all
% volumes

% startIdx : the index of the the run. If length(pointStats)=N, the total
% number of runs is N*nGroups. For a complete analysis, startIdx will take
% a value 1:N*nGroups for each run

% nRef : references to use for the recording. Will use up to nRef evenly
% spaced references.
% offset : added to start idx due to the limit on slurm jobs

%% initial parameters

%parameters for tracking/matching
param.dim=3;
param.good=2;
param.excessive=4;
param.quiet=1;
param.timeLimit=10;
param.difficult=1.5e4;

show=00;
%% default inputs
if nargin==0
    %if no input, manually select pointstats file
    filePath=uipickfiles;
    filePath=filePath{1};
    startIdx=1;
    start_offset=0;
elseif nargin==3
    start_offset=0;
end

%load pointStats file
load(filePath);
%% initial setup of which frames to select as reference and which to analyze
startIdx=startIdx+start_offset;
%each volume is analyzed ngroups times with each run having 150 references.
%The total number of references is ngroups*150.

%list of stacks presents
volList=find(cellfun(@(x) ~isempty(x),{pointStats.stackIdx}));
nVol=length(runIdxListAll);
presentIdx=[pointStats.stackIdx];


%% make output folders
outputFolder=fileparts(filePath);
outputFolder=[outputFolder filesep 'TrackMatrix'];

if ~isdir(outputFolder)
    mkdir(outputFolder)
end
%%




%the list of references to use
refList=unique(round(1:nRef/nVol:nVol));
refList=volList(refList);

%get sample points being matched
i_ps=presentIdx(startIdx);
P1=pointStats(i_ps);

%output path
outputName=[outputFolder filesep...
    'trackMatrix' num2str(i_ps,'%3.5d')...
    'Run' num2str(0,'%3.2d')];
display(outputName);

length_P1=size(P1.straightPoints,1);
%initialize trackMatrix, which will hold all the matches
TrackMatrixi=zeros(size(P1.straightPoints,1),nRef);
DMatrixi_x=TrackMatrixi;
DMatrixi_y=TrackMatrixi;
DMatrixi_z=TrackMatrixi;
%%
for iRef=1:nRef
    itic=tic;
    % get reference points being matched
    j_ps=refList(iRef);
    P2=pointStats(j_ps);
    try
        %%
        % get coordinates, volumes, and brightness for registration
        T1=[P1.straightPoints P1.Volume.^(1/3) P1.Rintensities];
        T2=[P2.straightPoints P2.Volume.^(1/3) P2.Rintensities];
        %make note of old points
        T1temp=P1.straightPoints(:,1:3);
        T2temp=P2.straightPoints(:,1:3);
        length_P2=size(P2.straightPoints,1);
        % do non rigid pointset regstration, first with entire pointset
        [T2_trans, ~, ~] = ...
            gmmreg_L2_multilevel_jn(...
            T2,T1,1, [ 1,.5], ...
            [.000008, 0.0000008, 0.0008],[0 0],...
            [ 0.000001 0.0001 0.001 0.001],show);
        
        %put together old points and transformed points for use in
        %kmeans clustering to break the pointsets into smaller groups
        
        trackInput=[T1temp ;T2_trans(:,1:3) ];
        idx = kmeans(trackInput(:,1:3),3);
        
        %get the cluster for each point in the sample and the reference
        idx1=idx(1:length_P1);
        idx2=idx(length_P1+1:end);
        %%
        % do non rigid pointset regstration with each of the cluseters
        % seperately
        
        for iRegions=1:max(idx)
            [T2_trans(idx2==iRegions,:), ~, ~] = ...
                gmmreg_L2_multilevel_jn(...
                T2_trans(idx2==iRegions,:),T1(idx1==iRegions,:), ...
                2, [ 3,.3,.3], ...
                [0.005,.0005, 0.002, 0.08],[0 0],...
                [ 0.000001 0.000001 0.000001 0.001],show);
        end
        
        %%
        track1=[];
        track2=[];
        %build track input, which contains:
        %1. the points to be matched (after transforming sample)
        %2. the untransformed points
        %3. the index of the point (1:length of points)
        %4. the time, either 1 or 2.
        trackInput_t1=[...
            T1temp T1temp  (1:length_P1)'  ones(length_P1,1)];
        trackInput_t2=[...
            T2_trans(:,1:3) T2temp  (1:length_P2)' 2*ones(length_P2,1)];
        
        trackInput=[trackInput_t1;trackInput_t2];
        %loop through each cluster and do tracking
        for iRegions=1:max(idx)
            %select only points in that cluster
            trackInputi=trackInput(idx==iRegions,:);
            %get different times of points, to make sure that both
            %time points are represented
            track_times=trackInputi(:,end);
            if length(unique(track_times))==2
                %counter is the max distance between matched points. If
                %the tracking fails, reduce this and try again.
                counter=18;
                TrackOut=nan;
                %do matching of points
                while(all(isnan(TrackOut(:))))
                    TrackOut=trackJN(trackInputi,counter,param);
                    counter=counter-1;
                end
                %% get the indices of the matched points
                TrackStats=round(TrackOut(:,7:end));
                track1i=TrackStats(1:2:end,1);
                track2i=TrackStats(2:2:end,1);
                track1=[track1;track1i];
                track2=[track2;track2i];
            end
        end
        %build track matrix, which shows at a given time, which points
        %were matched.
        TrackMatrixi(track1,runIdx)=track2;
        
        %%%%%%
        %playing around with distances, not currently used
        presentIJ=TrackMatrixi(:,runIdx)>0;
        points1=T1temp(track1,1:3);
        points2=T2_trans(track2,1:3);
        pointsDiff=abs(points1-points2);
        
        DMatrixi_x(presentIJ,runIdx)=pointsDiff(:,1);
        DMatrixi_y(presentIJ,runIdx)=pointsDiff(:,2);
        DMatrixi_z(presentIJ,runIdx)=pointsDiff(:,3);
        %%%%%%
        if ~mod(j_ps,10)
            disp(['Finished match' num2str(j_ps) ' in ' num2str((toc(itic))) 's']);
        end
    catch me
        display(me.identifier)
    end
end
if isempty(TrackMatrixi)
    TrackMatrixi=[];
end
%save TrackMatrixi
save(outputName,'TrackMatrixi');


