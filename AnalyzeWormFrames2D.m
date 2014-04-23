%function AnalyzeWormFrames2D(trackFolder)

%if nargin==0
    trackFolder=uigetdir;
%end
    maxDist=20;
    params.mem=4;
    params.dim=2;
    minTrack=1000;

files=dir([trackFolder filesep '*.mat']);
trackData=[];
trackIdx=0;
for imat=1:2:length(files)
    trackIdx=trackIdx+1;
    load([trackFolder filesep files(imat).name]);
    tracks=[centroids,Gintensities,Rintensities,Volume,trackIdx*ones(size(Volume))];
    trackData=[trackData;tracks];
    progressbar((1+imat)/length(files));
end

trackOutput=track(trackData,maxDist,params);
    trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
    
    badtracks=find(trackLengths<minTrack);
    badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
    
    trackOutput(badtracks,:)=[];
      %  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
         [ trackIdx,ia,ib]=unique(trackOutput(:,end));
         %%
for i=1:max(unique(ib))
             plot(trackOutput(ib==i,end-1),normalizeRange(smooth(trackOutput(ib==i,3),20)))%./ ...
                % normalizeRange(smooth(trackOutput(ib==i,4),100)));
             % plot(normalizeRange(smooth(trackOutput(ib==i,4)-mean(trackOutput(ib==i,4)),100)),'r');
             pause(1);
     %        hold off;
end


    

         

      

    
    
    
