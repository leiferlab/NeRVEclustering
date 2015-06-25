 dataFolder=uipickfiles;
 dataFolder=dataFolder{1};
 
 %%
rows=1200; cols=600;
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);


    vidInfo.bfAll=bfAll;
    vidInfo.fluorAll=fluorAll;
    vidInfo.hiResData=hiResData;
    

%%
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear','extrap');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));



stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
%stack2BFidx=stack2BFidx(~isnan(stack2BFidx))
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');


%%
zWave=hiResData.Z;
zWave=gradient(zWave);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,20));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));


%%

userList={'Fred' ,'Ashley','David','Mochi','Sagar','George','Jeff'};
NSegments=15;
stackLength=max(hiResData.stackIdx);
fiducialPointsTempAll=cell(stackLength,NSegments);
fiducialPoints=cell(stackLength,1);
pointStatsLookup={pointStats.stackIdx};
pointStatsLookup(cellfun(@(x) isempty(x), pointStatsLookup))={0};
pointStatsLookup=cell2mat(pointStatsLookup);
for i=1:length(pointStats)
    
        currentData=pointStats2(i);
if ~isempty(currentData.stackIdx)
    iStack=currentData.stackIdx;
        hiResIdx=find(hiResData.stackIdx==iStack);
    zRange=hiResData.Z(hiResIdx);
    
    rawPoints=currentData.rawPoints;
    trackIdx=currentData.trackIdx;
    if isfield(currentData,'trackWeights')
    trackWeights=currentData.trackWeights;
    end
    goodPoints=~isnan(trackIdx);
    rawPoints=rawPoints(goodPoints,:);
    trackIdx=trackIdx(goodPoints);
    trackWeights=trackWeights(goodPoints);
    zLevels=rawPoints(:,3);
   hiResIdx=hiResIdx+zOffset;
if 0%sign(nanmean(diff(zRange)))==-1
    
    zRange2=zRange([true ; diff(zRange)<0]);
    zInterp=interp1(flipud(hiResIdx),zLevels,'*linear');
zVoltages=interp1(flipud(zRange),zLevels,'*linear');
else
        zRange2=zRange([true ; diff(zRange)>0]);
    zInterp=interp1((hiResIdx),zLevels,'*linear');
zVoltages=interp1((zRange),zLevels,'*linear');

end

fiducialPointsTemp=cell(max(trackIdx),5);

fiducialPointsTemp(trackIdx,1:5)=num2cell([rawPoints(:,1:2),zVoltages,zInterp trackWeights]);
for iBot=1:NSegments
    fiducialPointsTempAll{iStack,iBot}=fiducialPointsTemp(iBot:NSegments:size(fiducialPointsTemp,1),:);
end

fiducialPoints{iStack-1}=fiducialPointsTemp;
end
    
end
fiducialPointsBot=fiducialPoints;
folderStr=[dataFolder filesep datestr(now,'yyyymmddTHHMMSS') 'Fiducials'];
timeOffset=zOffset;

%%
mkdir(folderStr)
clickPoints=0;
%fiducialPoints=fiducialPointsBot';
for iBot=1:NSegments
fiducialPoints=fiducialPointsTempAll(:,iBot);
fiducialPoints(cellfun(@(x) isempty(x),fiducialPoints))={cell(5,5)};
save([folderStr filesep 'Bot' num2str(iBot,'%3.2d')], 'fiducialPoints','clickPoints');

end
save([folderStr filesep 'timeOffset'], 'timeOffset');
% fiducialPoints=repmat({cell(5,5)},length(fiducialPointsBot),1);
% save([folderStr filesep 'Fred'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Ashley'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'David'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Mochi'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Sagar'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Andy'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'George'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Milos'], 'fiducialPoints','clickPoints');
% save([folderStr filesep 'Jeff'], 'fiducialPoints','clickPoints');
% 
% 
% 

