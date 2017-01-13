function [TrackMatrixAll,pointStats]=compileTrackMatrix(fileName)

%function takes the cluster output of multiple track matrix .mat files and
%compiles them into a cell array for further analysis.
if nargin==0
    fileName=uipickfiles;
    fileName=fileName{1};
end

dataFolder=fileName;
%make file list
fileList=dir([fileName filesep 'trackMatrix*.mat']);
if isempty(fileList);
    fileName=[fileName filesep 'TrackMatrix'];
    fileList=dir([ fileName filesep 'trackMatrix*.mat']);
end

%load pointstats file
pointStats=load([dataFolder filesep 'PointsStats.mat']);
pointStats=pointStats.pointStats;
%initialize TrackMatrixi field
for i=1:length(pointStats);
    pointStats(i).TrackMatrixi=[];
end

presentIdx=find(cellfun(@(x) any(x),{pointStats.stackIdx}));
if length(presentIdx)~=str2double(fileList(end).name(12:16));
    presentIdx=1:length(pointStats);
end
%% loop over files to load
for i=1:length(fileList)
    %naming structure hard coded
    currentFile=[fileName filesep fileList(i).name];
    display(['Loading file ' fileList(i).name]);
    currentFileIdx=str2double(fileList(i).name(12:16));
    currentPSIdx=presentIdx(currentFileIdx);
    if strfind(fileList(i).name,'Run')
        runIdx=1+str2double(fileList(i).name(20:21));
    else
        runIdx=1;
    end
    currentData=load(currentFile);
    %populate fields from individual pointStatsfiles
    for jname=fieldnames(currentData)'
        %    TrackMatrixAll{presentIdx(i)}=currentData.TrackMatrixi;
        if runIdx==1
            pointStats(currentPSIdx).(jname{1})=currentData.(jname{1});
        else
            pointStats(currentPSIdx).(jname{1})=[pointStats(currentPSIdx).(jname{1}) ...
                currentData.(jname{1})];
        end
        
    end
end
%no longer used
TrackMatrixAll=[];