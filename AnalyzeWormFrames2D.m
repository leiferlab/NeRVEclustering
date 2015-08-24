%function AnalyzeWormFrames2D(trackFolder)

%if nargin==0
    trackFolder=uigetdir;
%end
    maxDist=3;
    params.mem=40;
    params.dim=2;
    minTrack=1000;
    params.good=2;
    params.quiet=1;
files=dir([trackFolder filesep 'output*.mat']);
nfiles=(length(files));
% strLength=cellfun(@(x) length(x), {files.name});
% files=files(strLength==mode(strLength));

trackData=[];
trackIdx=0;
%%
close all
for imat=round(1:30:length(files))
    fileName=['output' num2str(imat,'%5.6d') '.mat'];
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
 %   pause(2);
end

%%
matSearch=round(1:nfiles);
%%
tracksi=repmat({[]},1,length(matSearch));

parfor i=1:length(matSearch)
    imat=matSearch(i);
  fileName=[trackFolder filesep 'output' num2str(imat,'%5.6d') '.mat'];
  if exist(fileName,'file')
      data=load(fileName);
     data=data.outputDataTemp;
    centroids=data.centroids;
    Gintensities=data.Gintensities-90;
    Rintensities=data.Rintensities-90;
    Volume=data.Volume;
    trackIData=[centroids,Gintensities,Rintensities,Volume,i*ones(size(Volume))];
 %   trackIData=trackIData(Rintensities>50,:);
    tracksi{i}=trackIData;
    
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
%refIdx=find(nPoints==nRef);
refIdx=100;
refData=tracksi{refIdx(end)};
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
%%
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

end
end
if ~mod(i,100)
display(['Finished ' num2str(i)]);
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
   % good=t>10;
    %t=t(good);a=a(good);b=b(good);
    a=normalizeRange(a);
    b=normalizeRange(b);
%    t=t/13;
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
%%

    amat(mean(~isnan(bmat),2)<.75,:)=[];
    bmat(mean(~isnan(bmat),2)<.75,:)=[];

         GvalsAll=amat;
RvalsAll=bmat;

      
%%
    
    
    
    photoBleachingR=zeros(size(RvalsAll));
photoBleachingG=zeros(size(RvalsAll));

%Fexponent =@(x,xdata) x(1)*exp(x(2)*xdata)+x(3);
Fexponent=fittype('a*exp(b*x)+c','dependent',{'y'},'independent',...
    {'x'},'coefficients',{'a', 'b', 'c'});
fitOptions=fitoptions(Fexponent);
fitOptions.Lower=[0,-.2,0];
fitOptions.Upper=[1001,0,10000];

progressbar(0)


for i=1:size(RvalsAll,1)
    
    try
        %%
     
    progressbar(i/size(RvalsAll,1));
    xVals=round(linspace(1,size(RvalsAll,2),10000));
    
    present=(~isnan(RvalsAll(i,xVals)+GvalsAll(i,xVals))') ;
    xVals=xVals(present);
    rVals=RvalsAll(i,xVals)';
    minWindow=2500;
    gVals=GvalsAll(i,xVals)';
    gVals=ordfilt2(gVals,1,true(minWindow,1));
    rVals=ordfilt2(rVals,1,true(minWindow,1));
    gVals(1:(minWindow/2))=gVals((minWindow/2));
    rVals(1:(minWindow/2))=rVals((minWindow/2));
    gVals(end-(minWindow/2):end)=gVals(end-(minWindow/2));
    rVals(end-(minWindow/2):end)=rVals(end-(minWindow/2));
       fitOptions.StartPoint=[1,-.1,min(rVals)];


f=fit(xVals,rVals,Fexponent,fitOptions);
   fitOptions.StartPoint=[1,-.1,min(gVals)];

g=fit(xVals,gVals,Fexponent,fitOptions);
if f(1)>(max(RvalsAll(i,:))+100)
    f=fit(xVals,rVals,'poly1');
if f.p1>0
    f.p1=0;
end
    
end
if g(1)>(max(GvalsAll(i,:))+100)
    g=fit(xVals,gVals,'poly1');
if g.p1>0
    g.p1=0;
end

  %  g.a=0;
end
%[f.b g.b]
% plot(GvalsAll(i,:))
% hold on
% plot(g)
% imagesc([],[0 1000],ethoTrack')
% alpha(.2);
% ylim([0 max(GvalsAll(i,:))]);
% hold off
% subplot(2,1,2);
% plot(RvalsAll(i,:))
% hold on
% plot(f)
% imagesc([],[0 1000],ethoTrack')
% alpha(.2);
% ylim([0 max(RvalsAll(i,:))]);
% hold off
% drawnow

photoBleachingR(i,:)=f((1:size(RvalsAll,2)))-f(size(RvalsAll,2));

photoBleachingG(i,:)=g((1:size(RvalsAll,2)))-g(size(RvalsAll,2));
    catch
    end
    

end
% photoBleachingR=bsxfun(@minus,photoBleachingR,photoBleachingR(:,end));
% photoBleachingG=bsxfun(@minus,photoBleachingG,photoBleachingG(:,end));
Rvalstemp=RvalsAll-photoBleachingR ;
Rvalstemp(bsxfun(@le,Rvalstemp,quantile(Rvalstemp,.1,2)))=nan;
Rvalstemp=colNanFill(Rvalstemp')';

Gvalstemp=GvalsAll-photoBleachingG ;
Gvalstemp(bsxfun(@le,Gvalstemp,quantile(Gvalstemp,.1,2)))=nan;
Gvalstemp=colNanFill(Gvalstemp')';

%%

RatioVals=Gvalstemp./Rvalstemp;
RvalsPhotoCorr=Rvalstemp;
GvalsPhotoCorr=Gvalstemp;

dataFolder=fileparts(trackFolder);


%%
save([dataFolder filesep 'singlePlaneData'],'RatioVals','RvalsPhotoCorr','GvalsPhotoCorr',...
    'RvalsAll','GvalsAll');
    
    
