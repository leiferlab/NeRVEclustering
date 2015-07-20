%function AnalyzeWormFrames2D(trackFolder)

%if nargin==0
    trackFolder=uigetdir;
%end
    maxDist=3;
    params.mem=40;
    params.dim=2;
    minTrack=1000;
    params.good=2;
files=dir([trackFolder filesep 'output*.mat']);
nfiles=min(length(files), 99999);
% strLength=cellfun(@(x) length(x), {files.name});
% files=files(strLength==mode(strLength));

trackData=[];
trackIdx=0;
%%
close all
for imat=round(1:100:length(files))
    fileName=['output' num2str(imat,'%5.5d') '.mat'];
      data=load([trackFolder filesep fileName]);
    data=data.outputDataTemp;
    centroids=data.centroids;
   % centroids=bsxfun(@minus,centroids,mean(centroids,1));
    if imat==1
    h=scatter(centroids(:,1),centroids(:,2));
         xlim(minmax(centroids(:,1)'));
     ylim(minmax(centroids(:,2)'));axis equal;
     else
         h.XData=centroids(:,1);
         h.YData=centroids(:,2);
     end
    drawnow
    pause(2);
end

%%
matSearch=round(1:nfiles);
%%
tracksi=repmat({[]},1,length(matSearch));

parfor i=1:length(matSearch)
    imat=matSearch(i);
  fileName=[trackFolder filesep 'output' num2str(imat,'%5.5d') '.mat'];
  if exist(fileName,'file')
      data=load(fileName);
     data=data.outputDataTemp;
    centroids=data.centroids;
    Gintensities=data.Gintensities;
    Rintensities=data.Rintensities;
    Volume=data.Volume;
    tracksi{i}=[centroids,Gintensities,Rintensities,Volume,i*ones(size(Volume))];
   % progressbar((1+imat)/length(files));
   display(['completed frame' num2str(i)]);
  end
    
end
%%
trackData=cell2mat({tracksi{:}}');

%%
nTimes=length(tracksi);
nPoints=cellfun(@(x) length(x), tracksi);
nRef=max(nPoints);
refIdx=10000;%find(nPoints==nRef);
refData=tracksi{refIdx};
nRef=size(refData,1);
refData(:,end)=1;
refData(:,5)=1:length(refData);
refCentroids=refData(:,1:2);
tracksOutputi=repmat({[]},1,nTimes);
offsetAll=[0 0];
show=1;
if ~show
    progressbar(0)

end

offsetAll=[-1 10];
for i=1:nTimes;
currentData=tracksi{i};
if ~isempty(currentData)
currentData(:,end)=2;
currentData(:,5)=1:length(currentData);
currentData(:,1:2)=bsxfun(@minus,currentData(:,1:2),offsetAll);

currentCentroids=currentData(:,1:2);
%indexPairs = matchFeatures(refCentroids,currentCentroids,'Unique',1,'MaxRatio',1,'MatchThreshold',100);
% offset=currentCentroids(indexPairs(:,2),:)-refCentroids(indexPairs(:,1),:);
% offset=median(offset,1);
% offsetAll=offset+offsetAll;

TrackStats=trackJN([refData;currentData],10,params);


TrackedIDs=TrackStats([1;diff(TrackStats(:,end))]==0,end);
 TrackStats=TrackStats(ismember(TrackStats(:,end),TrackedIDs),:);
track1=TrackStats(1:2:end,end-2);
track2=TrackStats(2:2:end,end-2);
indexPairs=[track1 track2];
offset=currentCentroids(track2,:)-refCentroids(track1,:);
offset=median(offset,1);
offsetAll=offset+offsetAll;


currentData(:,end)=i;
outputData=[currentData(indexPairs(:,2),:), indexPairs(:,1)];
tracksOutputi{i}=outputData;

if show
if ~mod(i,100)
          
    h=scatter(currentCentroids(:,1),currentCentroids(:,2));
    hold on
    scatter(refData(:,1),refData(:,2),'rx');
    htext=text(currentCentroids(indexPairs(:,2),1),...
        currentCentroids(indexPairs(:,2),2),...
        cellstr(num2str(indexPairs(:,1)))');
    
    htext2=text(refData(:,1)-10,...
        refData(:,2),...
        cellstr(num2str((1:nRef)')),'color','r');
    
     axis equal;
        xlim([50 500]);ylim([50 500]);
     hold off
          title(['Frame ' num2str(i)]);

     drawnow
%      else
%     h=scatter(currentCentroids(:,1),currentCentroids(:,2));
%     htext=text(currentCentroids(indexPairs(:,2),1),...
%         currentCentroids(indexPairs(:,2),2),...
%         cellstr(num2str(indexPairs(:,1)))');
%          xlim(minmax(currentCentroids(:,1)'));
%      ylim(minmax(currentCentroids(:,2)'));axis equal;
end
else
        progressbar(i/nTimes)

end
end
end
trackOutput=double(cell2mat(tracksOutputi'));
display('done!');

%%
% trackOutput=[];
% tryCounter=1;
% while isempty(trackOutput)
% trackOutput=trackJN(trackData,maxDist/tryCounter,params);
% tryCounter=tryCounter+1;
% 
% end
% 
%     trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
%     
%     badtracks=find(trackLengths<minTrack);
%     badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
%     
%     trackOutput(badtracks,:)=[];
      %  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
         %%
                  [ trackIdx,ia,ib]=unique(trackOutput(:,end));
show=1;
amat=nan(max(ib),max(trackOutput(:,end-1)));
bmat=amat;
for i=1:max(unique(ib))
    t=trackOutput(ib==i,end-1);
    a=((trackOutput(ib==i,3)));
    b=((trackOutput(ib==i,4)));
    amat(i,t)=a;
    bmat(i,t)=b;
    good=t>10;
    t=t(good);a=a(good);
    a=normalizeRange(a);
    b=normalizeRange(b(good));
    t=t/13;
    if show
             plot(t,a,'g')%./ ...
             hold on
             plot(t,b,'r')
             hold off
                % normalizeRange(smooth(trackOutput(ib==i,4),100)));
             % plot(normalizeRange(smooth(trackOutput(ib==i,4)-mean(trackOutput(ib==i,4)),100)),'r');
             pause(.1);
     %        hold off;
    end
end


    

         

      

    
    
    
