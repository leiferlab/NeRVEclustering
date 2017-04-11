function clusterWormTracker(dataFolder,startIdx,stepSize)
% clusterWormTracker compares a set of pointsets from WormStraighten code
% and uses non rigid pointset registration to create match matrices. The
% code has been now modified to just run with on one sample with nRef for a
% given pointStats. You can no longer split up a single sample.

%%%% Inputs
% filePath : link to the complete PointStats file with data from all
% volumes, also needs to have a PointStatsRef file which has the data from
% the volumes selected to be references

% startIdx : the index of the volume being analyzed. Only volumes that are
% not empty are analyzed. 



%% default inputs
if nargin==0
    %if no input, manually select pointstats file
    dataFolder=uipickfiles;
    dataFolder=dataFolder{1};
    startIdx=1;
end

if nargin==2
    stepSize=1;
end


PS_file=[dataFolder filesep 'PointsStats.mat'];
ref_file=[dataFolder filesep 'pointStatsRef.mat'];
%load pointStats file
load(PS_file);
%load reference file
load(ref_file)


%% make output folders
outputFolder=[dataFolder filesep 'TrackMatrix'];

if ~isdir(outputFolder)
    mkdir(outputFolder)
end
%% initial setup of which frames to select as reference and which to analyze
%each volume is analyzed ngroups times with each run having 150 references.
%The total number of references is ngroups*150.


presentIdx=[pointStats.stackIdx];

for iRun=1:stepSize
    %get sample points being matched
    currentIdx=startIdx+iRun-1;
    i_ps=presentIdx(currentIdx);
    P1=pointStats(i_ps);
    
    %output path
    outputName=[outputFolder filesep...
        'trackMatrix' num2str(i_ps,'%3.5d')...
        'Run' num2str(0,'%3.2d')];
    display(outputName);
    
    
    %%
    % do matching with references
    TrackMatrixi=compareWithRef(P1,PS_ref);
    if isempty(TrackMatrixi)
        TrackMatrixi=[];
    end
    %save TrackMatrixi
    save(outputName,'TrackMatrixi');
end


