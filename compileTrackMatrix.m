function [TrackMatrixAll,pointStats]=compileTrackMatrix(fileName)

%function takes the cluster output of multiple track matrix .mat files and
%compiles them into a cell array for further analysis.
if nargin==0
    fileName=uipickfiles;
    fileName=fileName{1};
end

dataFolder=fileName;

fileList=dir([fileName filesep 'trackMatrix*.mat']);
if isempty(fileList);
    fileName=[fileName filesep 'trackMatrix'];
    fileList=dir([ fileName filesep 'trackMatrix*.mat']);
end


pointStats=load([dataFolder filesep 'PointsStats']);
pointStats=pointStats.pointStats;

for i=1:length(pointStats);
    pointStats(i).TrackMatrixi=[];
end

for i=1:length(fileList)
    currentFile=[fileName filesep fileList(i).name];
    display(['Loading file ' fileList(i).name]);
    currentFileIdx=str2double(fileList(i).name(12:16));
    if strfind(fileList(i).name,'Run')
        runIdx=1+str2double(fileList(i).name(20:21));
    else
        runIdx=1;
    end
    currentData=load(currentFile);
    for jname=fieldnames(currentData)'
        %    TrackMatrixAll{presentIdx(i)}=currentData.TrackMatrixi;
        if runIdx==1
            pointStats(currentFileIdx).(jname{1})=currentData.(jname{1});
        else
            pointStats(currentFileIdx).(jname{1})=[pointStats(currentFileIdx).(jname{1}) ...
                currentData.(jname{1})];
        end
        
    end
end
TrackMatrixAll=[];