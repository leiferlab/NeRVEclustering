dataFolder=uipickfiles;
dataFolder=dataFolder{1};
imSize=[1200 600];
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


%%
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');

%%
behavior=load([dataFolder filesep 'ManualBehavior']);
behavior=behavior.behavior;

hiResBehavior=behavior(diff(BF2stackIdx)>0);

%%

xPos=hiResData.xPos(diff(hiResData.stackIdx)==1);
yPos=hiResData.yPos(diff(hiResData.stackIdx)==1);
xPos=xPos(hasPoints);
yPos=yPos(hasPoints);


