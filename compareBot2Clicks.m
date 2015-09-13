%load pointStatsFiles and coordinate transformation lookup tables to check
%accuracy of automated tracking. coordinate transforms are currently saved
%in individual pointstats files

imageFolder='V:\20150616\BrainScanner20150616_100353-worm1\CLstraight_20150813\';
presentIdx=1:length(pointStats2);
param.dim=3;
param.excessive=4;
 param.quiet=1;
 param.difficult=2.e4;
 param.good=2;
 %%
for i=1:length(pointStats2)
   %      pointStats2(presentIdx(i)).transitionMatrix=TrackMatrix{i};
   
    try
        %%
        iFile=pointStats2(i).stackIdx;
           pointFile=([imageFolder  filesep 'pointStats' num2str(iFile,'%2.5d') '.mat']);

     pointStatsTemp=load(pointFile);
pointStatsTemp=pointStatsTemp.pointStats;
X=double(pointStatsTemp.transformx);
Y=double(pointStatsTemp.transformy);
Z=double(pointStatsTemp.transformz);

    rawPoints=pointStats2(presentIdx(i)).straightPoints;
    rawPoints=rawPoints(~isnan(pointStats2(presentIdx(i)).trackIdx),:);
    rawPointsId=pointStats2(presentIdx(i)).trackIdx(~isnan(pointStats2(presentIdx(i)).trackIdx));
    rawPointsId=(1:length(rawPointsId))';
    

    rawPoints=coordinateTransform3d(rawPoints,X,Y,Z);
%         present=~isnan(rawPoints(:,1));
%     present=find(present);
%     rawPointsId=rawPointsId(present);
%     rawPoints=rawPoints(present,:);
%     
    
    controlPoints=pointStats2(presentIdx(i)).controlPoints;
    rawPoints=[rawPoints rawPointsId ones(size(rawPoints,1),1)];
    
    controlPoints=[controlPoints (1:size(controlPoints,1))' 2*ones(size(controlPoints,1),1)];
    matchOut=trackJN([rawPoints;controlPoints],9,param);
    matchID=matchOut(:,end-2);
    pairedId=reshape(matchOut(:,end-2),2,[])';
    trackIdx=pointStats2(presentIdx(i)).trackIdx;
    trackIdx(~ismember(trackIdx,pairedId(:,1)))=nan;
     matchIdx=nan(size(pointStats2(presentIdx(i)).trackIdx));
    matchIdx((pairedId(:,1)))=pairedId(:,2);
    
    pointStats2(presentIdx(i)).matchIdx=matchIdx;
    display(['Finished frame : ' num2str(presentIdx(i))]);
    catch ME
        display(['Error frame:' num2str(presentIdx(i))]);
    end
    
end
