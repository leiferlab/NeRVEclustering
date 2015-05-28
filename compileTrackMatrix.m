function [TrackMatrixAll,pointStats]=compileTrackMatrix(fileName)

%function takes the cluster output of multiple track matrix .mat files and
%compiles them into a cell array for further analysis.

if nargin==0
    fileName=uipickfiles;
    fileName=fileName{1};
end



fileList=dir([fileName filesep 'trackMatrix*.mat']);


pointStats=load([fileName filesep 'PointsStats']);
pointStats=pointStats.pointStats;
% presentIdx=cellfun(@(x) ~isempty(x),{pointStats.stackIdx},'uniform',0);
% presentIdx=find(cell2mat(presentIdx));


for i=1:length(fileList)
    currentFile=[fileName filesep fileList(i).name];
    currentFileIdx=str2double(currentFile(end-8:end-4));
    currentData=load(currentFile);
    TrackMatrixAll{i}=currentData.TrackMatrixi;
    pointStats(i).TrackMatrixi=currentData.TrackMatrixi;
    
end
