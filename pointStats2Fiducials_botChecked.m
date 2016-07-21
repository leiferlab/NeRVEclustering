pointStatsFile=uipickfiles();
pointStatsFile=pointStatsFile{1};
dataFolder=fileparts(pointStatsFile);

imageFolder=uipickfiles('filterspec',dataFolder);
imageFolder=imageFolder{1};
load(pointStatsFile)

%%
pointStats2=pointStatsNew;

presentIdx=1:length(pointStats2);
rows=1200;
cols=600;
%%
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);

zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
timeOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%%

fiducialPoints=cell(1,length(pointStats2));
for i=1:length(pointStats2)
   %      pointStats2(presentIdx(i)).transitionMatrix=TrackMatrix{i};
   
    try
        %%
        iFile=pointStats2(i).stackIdx;
        if ~isempty(iFile)
        pointFile=([imageFolder  filesep 'pointStats' num2str(iFile,'%2.5d') '.mat']);

     pointStatsTemp=load(pointFile);
pointStatsTemp=pointStatsTemp.pointStats;
X=double(pointStatsTemp.transformx);
Y=double(pointStatsTemp.transformy);
Z=double(pointStatsTemp.transformz);
trackIdx=pointStats2(i).trackIdx;
    rawPoints=pointStats2(presentIdx(i)).straightPoints;
  %  rawPoints=rawPoints(~isnan(pointStats2(presentIdx(i)).trackIdx),:);
    present=~isnan(trackIdx+rawPoints(:,1));
    trackIdx=trackIdx(present);
    rawPoints=rawPoints(present,:);
    
    
    rawPointsId=pointStats2(presentIdx(i)).trackIdx(~isnan(pointStats2(presentIdx(i)).trackIdx));
    rawPointsId=(1:length(rawPointsId))';
    

    rawPoints2=coordinateTransform3d(rawPoints,X,Y,Z);
    rawPoints2(rawPoints2==0)=1;
    hiResRange=find(hiResData.stackIdx==pointStats2(i).stackIdx);
    hiResIdx=hiResRange(rawPoints2(:,3))+timeOffset;
    rawPoints2(rawPoints2==1)=nan;
    hiResIdx(isnan(rawPoints2(:,3)))=3;
    hiResV=interp1(hiResData.Z,hiResIdx-timeOffset);
    fiducialPointsi=cell(200,4);
    
    fiducialInput=num2cell([rawPoints2(:,1:2) hiResV hiResIdx]);
 %   rawPointsCell=num2cell(rawPoints);

     fiducialPointsi(trackIdx,:)=fiducialInput;
%     
%      present=~isnan(rawPoints(:,1));
%      present=find(present);
%      
%      rawPointsId=rawPointsId(present);
%      rawPoints=rawPoints(present,:);
% %     
%     
%     controlPoints=pointStats2(presentIdx(i)).controlPoints;
%     rawPoints=[rawPoints rawPointsId ones(size(rawPoints,1),1)];
%     controlPoints=[controlPoints (1:size(controlPoints,1))' 2*ones(size(controlPoints,1),1)];
%     matchOut=trackJN([rawPoints;controlPoints],9,param);
%     matchID=matchOut(:,end-2);
%     pairedId=reshape(matchOut(:,end-2),2,[])';
%     trackIdx=pointStats2(presentIdx(i)).trackIdx;
%     trackIdx(~ismember(trackIdx,pairedId(:,1)))=nan;
%      matchIdx=nan(size(pointStats2(presentIdx(i)).trackIdx));
%     matchIdx((pairedId(:,1)))=pairedId(:,2);
%     
%     pointStats2(presentIdx(i)).matchIdx=matchIdx;
     fiducialPoints{iFile}=fiducialPointsi;
    display(['Finished frame : ' num2str(presentIdx(i))]);
        end
    catch ME
        rethrow(ME)
    end
    
end


%%
newFiducialFolder=[dataFolder filesep 'BotfFiducialPoints'];
mkdir(newFiducialFolder);
clickPoints=0;
save([newFiducialFolder filesep 'timeOffset'],'timeOffset');
save([newFiducialFolder filesep 'botFiducials'],'fiducialPoints','clickPoints');

