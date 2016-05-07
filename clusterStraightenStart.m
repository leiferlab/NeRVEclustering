
function clusterStraightenStart(dataFolder)
rows=1200;
cols=600;
nPix=rows*cols;


    [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);

%bundle timing data

    vidInfo.bfAll=bfAll;
    vidInfo.fluorAll=fluorAll;
    vidInfo.hiResData=hiResData;
    

%% make time alignments
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear','extrap');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear','extrap');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));

stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
stack2fluoridx=fluorIdxLookup(diff(hiResData.stackIdx)==1);

topLimit=min(max(hiResData.stackIdx),length(stack2BFidx));
BF2stackIdx=interp1(stack2BFidx,1:topLimit,bfIdxList,'nearest');
fluor2stackIdx=interp1(stack2fluoridx,1:topLimit,fluorIdxList,'nearest');


%% load centerline
[centerline, CLoffset]=loadCLBehavior(dataFolder);

%%

minStart=max([min(BF2stackIdx(CLoffset+1)) min(fluor2stackIdx)])+1;

%% load alignment data

    %try to load the alignment file in the folder, otherwise, select them
    %individual in the registration folder
    alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;


%% calculate offset between frame position and zposition
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%% Create straightened destination folder
destination= ['CLstraight_' datestr(now,'yyyymmdd')];
imageFolder2=[dataFolder filesep destination];
mkdir(imageFolder2);

%% Select range of volumes to analyze

startStack=minStart;
endStack=max(hiResData.stackIdx);
stackRange= startStack:endStack;
%% do a sample image in range, for observation and for reference
tic
show=1;

counter=50; %which volume to do

%Run straighten and segmentation on one volume
[V,pointStatsOut,Vtemplate,side,lastOffset,Vbw]=...
    WormCLStraighten_10(dataFolder,destination,vidInfo,...
    alignments,[],[],zOffset,minStart+counter,'left',[],show);


% fields in the pointstats file that will be used
poinStatsFields={'straightPoints','rawPoints','stackIdx','pointIdx',...
    'Rintensities','Volume','controlPoints'};

for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=pointStatsOut.(field);
end

show=0;
display(['Finished image ' num2str(startStack,'%3.5d') ' in ' num2str(toc) 's'])
%save initial workspace for testing and later use
save([dataFolder filesep 'startWorkspace'])
