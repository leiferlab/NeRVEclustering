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
try
trackOutput=track(trackData,maxDist,params);
catch 
    trackOutput=track(trackData,maxDist/2,params);

end

    trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
    
    badtracks=find(trackLengths<minTrack);
    badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
    
    trackOutput(badtracks,:)=[];
      %  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
         [ trackIdx,ia,ib]=unique(trackOutput(:,end));
         %%
for i=1:max(unique(ib))
    t=trackOutput(ib==i,end-1);
    a=(smooth(trackOutput(ib==i,3),20));
    b=(smooth(trackOutput(ib==i,4),20));
    good=t>10;
    t=t(good);a=a(good);
    a=normalizeRange(a);
    b=normalizeRange(b(good));
    t=t/13;
             plot(t,a,'g')%./ ...
             hold on
             plot(t,b,'r')
             hold off
                % normalizeRange(smooth(trackOutput(ib==i,4),100)));
             % plot(normalizeRange(smooth(trackOutput(ib==i,4)-mean(trackOutput(ib==i,4)),100)),'r');
             pause(1);
     %        hold off;
end


    

         

      

    
    
    
