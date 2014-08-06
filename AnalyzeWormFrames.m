function AnalyzeWormFrames(trackFolder)
%%
if nargin==0
    trackFolder=uigetdir;
end
%%
    maxDist=20;
    params.mem=40;
params.good=100;
files=dir([trackFolder filesep '*.mat']);
trackData=[];
trackIdx=0;
Volume=[];
for imat=1:1:length(files)
    trackIdx=imat;
    load([trackFolder filesep files(imat).name]);
    centroids(:,3)=centroids(:,3)-mean(centroids(:,3));
%     scatter3(centroids(:,1),centroids(:,2),centroids(:,3),'.g')
%     hold on

tracks=[centroids,Gintensities,Rintensities,Volume,(1:length(Rintensities))',trackIdx*ones(size(Volume))];
trackTitles={'x','y','z','green','red','Volume','Frame Number', 'inFrame Index','trackIdx'}  ; 
trackData=[trackData;tracks];
   % progressbar((1+imat)/length(files));
end
%progressbar(1);
%%
    params.dim=size(centroids,2);
    success=0;
    maxDistFactor=1;
    while ~success
try
trackOutput=track(trackData,maxDist*maxDistFactor,params);
success=1;
catch
maxDistFactor=maxDistFactor*.7;
end
    end
    
    trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
    

      %  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
         [ trackIdx,ia,ib]=unique(trackOutput(:,end));
         
         save([trackFolder filesep 'trackOutput'],'trackOutput','trackTitles');
         %% plot intensities
for i=1:max(unique(ib))
           %  plot(trackOutput(ib==i,params.dim+1)./(trackOutput(ib==i,params.dim)),'g');
             plot(trackOutput(ib==i,params.dim+1),'g');
%              hold on
%              plot(),'r');
             pause(1);
             hold off;
end

         %% plot all positions of tracked particles
for i=1:max(unique(ib))
           %  plot(trackOutput(ib==i,params.dim+1)./(trackOutput(ib==i,params.dim)),'g');
             cline(trackOutput(ib==i,1),trackOutput(ib==i,2),trackOutput(ib==i,3),...
                 trackOutput(ib==i,end-1));
             axis equal
             hold on
             pause(.01);
%              plot(),'r');
end


%%
startTime=0;
normalizeFlag=0;
smoothWindow=5;
nTracks=max(trackOutput(:,end));
nTime=max(trackOutput(:,end-1));
activityMat=zeros(nTracks,nTime);
output=trackOutput(:,4);
for i=1:nTracks
    t=trackOutput((trackOutput(:,end)==i),end-1);
    a=output((trackOutput(:,end)==i));
    a=a(t>startTime);
    t=t(t>startTime);        
    a=(smooth(a,smoothWindow));
    if normalizeFlag
        a=normalizeRange(a);
    end
    
    activityMat(i,t)=a;
    
end
    
         
         
      

    
    
    
