function clusterWormStraightening(dataFolder,nStart,nRange)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calls worm straightening code if a startWorkspace.mat file is already
% created in the dataFolder being analyzed, program runs to straighten
% stacks nStart:nStart+nRange-1. Made to run on a cluster
% Inputs:
%   dataFolder- same requirements as the clusterStraightenStart program,
%               including the output of that program, startWorkspace.mat
%   nStart-     Volume to start straightening, made for cluster array jobs
%   nRange-     Number of volumes to straighten in this job.


%% load initial variables
straightenData=load([dataFolder filesep 'startWorkspace.mat']);

destination=straightenData.destination;
Vtemplate=straightenData.Vtemplate;
zOffset=straightenData.zOffset;
side=straightenData.side;
vidInfo=straightenData.vidInfo;

%% make output folders
display(dataFolder)
imageFolder2=[dataFolder filesep destination];

%% load alignments
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;

%% loop through range of stacks to do straightening
for iStack=nStart:(nStart+nRange-1)
    %set up image and pointstats names
    fileName2=[imageFolder2 filesep 'image' num2str(iStack,'%3.5d') '.tif'];
    fileName3=[imageFolder2 filesep 'pointStats' num2str(iStack,'%3.5d') '.mat'];
    % does not overwrite if both files are present
    tic
    if ~exist(fileName2,'file') && ~exist(fileName3,'file')
        % do straightening for this istack
        WormCLStraighten_11(dataFolder,destination,vidInfo,...
            alignments,Vtemplate,zOffset,iStack,side,0); 
        display(['image' num2str(iStack,'%3.5d') 'completed in ' num2str(toc) 's']);
    elseif exist(fileName3,'file')
        PS=load(fileName3);
        if isempty(PS.pointStats.straightPoints);
                 WormCLStraighten_11(dataFolder,destination,vidInfo,...
            alignments,Vtemplate,zOffset,iStack,side,0); 
        display(['image' num2str(iStack,'%3.5d') 'completed in ' num2str(toc) 's']);
        end
    else
        display([ 'image' num2str(iStack,'%3.5d') '.tif already exist!'])
    end
end
