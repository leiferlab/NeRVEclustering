dataFolder=uipickfiles;
dataFolder=dataFolder{1};
imSize=[1200 600];
%%

[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);


%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear','extrap');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


%%
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
%stack2BFidx=stack2BFidx(~isnan(stack2BFidx))
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');

%%
behavior=load([dataFolder filesep 'ManualBehavior']);
behavior=behavior.behavior;

hiResBehavior=behavior(diff(BF2stackIdx)>0);
area(hiResBehavior)

%%


