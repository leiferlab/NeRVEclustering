function [bfTime,hiResFrameTime,hasPoints,bfRange,hiResRange,hasPointsTime,lookup]=dataTimeAlignment(dataFolder,fiducialPoints)
%bfTime - 
imSize=[1200,600];
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime),'linear','extrap');
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)>=1);
nanList=isnan(stack2BFidx);
BF2stackIdx=interp1(stack2BFidx(~nanList),find(~nanList),bfIdxList,'nearest');
%%
lookup.BF2hiRes=bfIdxLookup;
lookup.Fluor2hiRes=fluorIdxLookup;
lookup.stack2BFidx=stack2BFidx;
lookup.BF2stackIdx=BF2stackIdx;

%%
% hasPoints=cellfun(@(x) size(cell2mat(x),1), fiducialPoints,'uniformoutput',0);
% hasPoints=cell2mat(hasPoints);
% hasPoints=find(hasPoints>max(hasPoints)*.9);
% hasPoints=hasPoints(hasPoints<max(BF2stackIdx));
 bfRange=[1 find(diff(BF2stackIdx)==1)];
% hasPoints=hasPoints(hasPoints<length(bfRange));
% hasPoints=min(hasPoints):max(hasPoints);
hasPoints=1:length(fiducialPoints);

hasPointsTime=hiResData.frameTime(diff(hiResData.stackIdx)==1);
hasPointsTime=hasPointsTime(hasPoints);
timeStart=min(hasPointsTime);
hasPointsTime=hasPointsTime-timeStart;
frameTime=hiResData.frameTime;

hiResRangeTrack=(find(diff(hiResData.stackIdx)>0));
hiResRangeTrack=hiResRangeTrack(hasPoints);
hiResRange=hiResRangeTrack(1):hiResRangeTrack(end);
hiResFrameTime=frameTime(hiResRange);
hiResFrameTime=hiResFrameTime-timeStart;


%bfRange=bfRange(hasPoints(1)):bfRange(hasPoints(end));
bfTime=bfAll.frameTime;
bfTime=bfTime-timeStart;
