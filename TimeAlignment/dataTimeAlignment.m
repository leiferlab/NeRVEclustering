function [bfTime,hasPoints,hasPointsTime]=dataTimeAlignment(dataFolder)
%bfTime - 
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);

%% load pointStats File
pointStatsFile=[dataFolder filesep 'pointStatsNew.mat'];
pointStats=load(pointStatsFile);

%%

hasPoints=1:length(pointStats);
hasPointsTime=hiResData.frameTime(diff(hiResData.stackIdx)==1);
hasPointsTime=hasPointsTime(hasPoints);

%bfRange=bfRange(hasPoints(1)):bfRange(hasPoints(end));
bfTime=bfAll.frameTime;
