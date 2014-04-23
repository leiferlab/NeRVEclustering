function AnalyzeWormFrames(trackFolder)

if nargin==0
    trackFolder=uigetdir;
end
    maxDist=20;
    params.mem=40;
    minTrack=100;

files=dir([trackFolder filesep '*.mat']);
trackData=[];
trackIdx=0;
for imat=1:2:8000%length(files)
    trackIdx=trackIdx+1;
    load([trackFolder filesep files(imat).name]);
    tracks=[centroids,Gintensities,Rintensities,Volume,trackIdx*ones(size(Volume))];
    trackData=[trackData;tracks];
    progressbar((1+imat)/length(files));
end
    params.dim=size(centroids,2);

trackOutput=track(trackData,maxDist,params);
    trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
    
    badtracks=find(trackLengths<minTrack);
    badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
    
    trackOutput(badtracks,:)=[];
      %  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
         [ trackIdx,ia,ib]=unique(trackOutput(:,end));
         %%
for i=1:max(unique(ib))
           %  plot(trackOutput(ib==i,params.dim+1)./(trackOutput(ib==i,params.dim)),'g');
             plot(trackOutput(ib==i,params.dim+2),'g');
%              hold on
%              plot(),'r');
             pause(1);
             hold off;
end


    
         
         
      

    
    
    
