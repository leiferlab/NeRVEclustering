function [TrackMatrixAll,pointStats]=compileTrackMatrix2(fileName)

%function takes the cluster output of multiple track matrix .mat files and
%compiles them into a cell array for further analysis.

if nargin==0
    fileName=uipickfiles;
    fileName=fileName{1};
end


fileList=dir([fileName filesep 'trackMatrix*Run*.mat']);
fileListCell={fileList.name};
trackStrIdx=(strfind(fileListCell{1},'x')+1):(strfind(fileListCell{1},'Run')-1);
runStrIdx=(strfind(fileListCell{1},'n')+1):(strfind(fileListCell{1},'.mat')-1);
trackIdxList=cellfun(@(x) str2double(x(trackStrIdx)),fileListCell);
runIdxList=cellfun(@(x) str2double(x(runStrIdx)),fileListCell);
[trackIdxListU,ia,ib]=unique(trackIdxList)

pointStats=load([fileName filesep 'PointsStats']);
pointStats=pointStats.pointStats;
 presentIdx=cellfun(@(x) ~isempty(x),{pointStats.stackIdx},'uniform',0);
 presentIdx=find(cell2mat(presentIdx));
presentIdx=1:max(presentIdx);
% presentIdx=cellfun(@(x) ~isempty(x),{pointStats.stackIdx},'uniform',0);
% presentIdx=find(cell2mat(presentIdx));


for i=trackIdxListU
    trackIdx=find(trackIdxList==i);
    currentTrackMatrix=0;
    for j=trackIdx
    currentFile=[fileName filesep fileList(j).name];
    currentFileIdx=j;
    currentData=load(currentFile);
    
   currentTrackMatrix=currentTrackMatrix+currentData.TrackMatrixi;
    end
    TrackMatrixAll{presentIdx(i)}=currentTrackMatrix;
    pointStats(presentIdx(i)).TrackMatrixi=currentTrackMatrix;
    
end
