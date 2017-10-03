function TrackMatrixi=compareWithRef(P1,PS_ref)


% CompareWithRef computes registration vectors between a single pointStats
% strucutre in P1 and an array of reference pointStats in PS_ref
% INPUTS:
%       P1 - pointStats structure from a single volume
%       PS_ref - pointStats array structure with a pointStats for all of
%       the reference volumes
% OUTPUTS:
%       TrackMatrix - A matrix of the all of the neuron registration
%       vectors from the neurons found in P1. TrackMatrix is an nxm matrix
%       where n is the number of neurons in P1 and m is the total number of
%       neurons in all of the volumes represented. 

%% initial parameters

%parameters for tracking/matching
param.dim=3;
param.good=2;
param.excessive=4;
param.quiet=1;
param.timeLimit=10;
param.difficult=1.5e4;
show=0;
nRef=length(PS_ref);
%%
TrackMatrixi=zeros(size(P1.straightPoints,1),nRef);
length_P1=size(P1.straightPoints,1);
%list of stacks presents

for j_ps=1:length(PS_ref)
    itic=tic;
    % get reference points being matched
    P2=PS_ref(j_ps);
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
        %This is to try and get a finer match using a smaller pointset. 
        
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
        TrackMatrixi(track1,j_ps)=track2;

        if ~mod(j_ps,10)
            disp(['Finished match' num2str(j_ps) ' in ' num2str((toc(itic))) 's']);
        end
    catch me
        display(me.identifier)
    end
end
