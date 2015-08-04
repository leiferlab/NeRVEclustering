pointStatsNew=pointStats2;
submissionFolder=uipickfiles('filterSpec','Y:\Jeff\');

submissionFolder=submissionFolder{1};
fileList=dir([submissionFolder filesep 'bot*.mat']);
pointStatsFile=dir([submissionFolder filesep 'point*']);
%%
load([submissionFolder filesep pointStatsFile.name]);
%%


    iPoint=str2double(fileList(1).name(11:15));
    data=load([submissionFolder filesep fileList(i).name]);
    
    xmean=nanmean(data.comparePointEstimate_x)';
    
  newXAll=nan(length(fileList),length(xmean));
  newYAll=newXAll;
  newZAll=newXAll;
  zScoreAll=newXAll;
  compareAllX=zeros(size(data.comparePointEstimate_x,1),...
      size(data.comparePointEstimate_x,2),length(fileList));
  compareAllY=compareAllX;
  compareAllZ=compareAllX;
  
parfor i=1:length(fileList)
    fileName=['botChecker' num2str(i,'%3.5d') '.mat'];
    display(['loading ' fileName]);
if exist([submissionFolder filesep fileName],'file')
    data=load([submissionFolder filesep fileList(i).name]);
    comparePointEstimate_x=data.comparePointEstimate_x;
    comparePointEstimate_y=data.comparePointEstimate_y;
    comparePointEstimate_z=data.comparePointEstimate_z;
    xmean=nanmean(comparePointEstimate_x)';
    xstd=nanstd(comparePointEstimate_x)';
    ymean=nanmean(comparePointEstimate_y)';
    ystd=nanstd(comparePointEstimate_y)';
    zmean=nanmean(comparePointEstimate_z)';
    zstd=nanstd(comparePointEstimate_z)';
    compareAllX(:,:,i)=comparePointEstimate_x;
    compareAllY(:,:,i)=comparePointEstimate_y;
    compareAllZ(:,:,i)=comparePointEstimate_z;
    xyzRefAll_zscore=data.xyzRefAll;
    xyzRefAll_zscore=xyzRefAll_zscore-[xmean ymean zmean];
    xyzRefAll_zscore=xyzRefAll_zscore./[xstd ystd zstd];
    zDistance=sqrt(sum(xyzRefAll_zscore.^2,2));
    
    check=find(zDistance>3);
    
    newX=nanmedian(comparePointEstimate_x);
    newY=nanmedian(comparePointEstimate_y);
    newZ=nanmedian(comparePointEstimate_z);
    
    newXAll(i,:)=newX;
    newYAll(i,:)=newY;
    newZAll(i,:)=newZ;
    
    zScoreAll(i,:)=zDistance;
end
end
[~,~,colorAll]=ndgrid(1:200,1:1522,1:93);


%%

check=(zScoreAll>4);

for i=1:length(pointStatsNew);
    if any(check(:,i))
        replaceIdx=find(check(:,i));
        
        for iReplace=replaceIdx'
            lookupIdx=pointStatsNew(i).trackIdx==iReplace;
            if any(lookupIdx)
            pointStatsNew(i).straightPoints(lookupIdx,:)=...
                [newXAll(iReplace,i),newYAll(iReplace,i),newZAll(iReplace,i)];
            else
                
            end
            
            end    
        
    end
end




%%

newXAll(newXAll==0)=nan;
newYAll(newXAll==0)=nan;
newZAll(newZAll==0)=nan;

for i=1:length(pointStatsNew);
    pointStatsNew(i).straightPoints=[newXAll(:,i),newYAll(:,i),newZAll(:,i)];
pointStatsNew(i).trackIdx=1:size(newXAll,1);
pointStatsNew(i).stackIdx=i;

end

