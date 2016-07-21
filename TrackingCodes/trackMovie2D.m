
    imFolder=uigetdir;
imFiles=dir([imFolder filesep  '*.tif']);
trackFolder=[imFolder filesep 'stackData'];
    maxDist=10;
    params.mem=10;
    minTrack=1000;

files=dir([trackFolder filesep '*.mat']);
load([imFolder filesep 'stackInfo.mat']);

%%
trackData=[];
trackIdx=0;
for imat=1:1:length(files)
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
         for iImage=1:10:length(files)
        wormMask=load([trackFolder filesep files(iImage).name],'wormMask');
wormMask=wormMask.wormMask;
        temp=double(imread([imFolder filesep imFiles(iImage).name],'tif'));
        temp=pixelIntensityCorrection(temp);
        temp_activity=temp((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
        worm=temp((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
        temp_activity=imwarp(temp_activity,t_concord,'OutputView',Rsegment);
        temp_activity(padRegion)=median(temp_activity(~padRegion));
        activity=bpass_jn(temp_activity,1,[40,40]);
        
         plotTrackData(activity,trackOutput,iImage);
         hold on
         B=bwboundaries(wormMask);
         for i=1:length(B)
             b=B{i};
             plot(b(:,2),b(:,1),'w')
         end
         hold off
         pause(.02)
         end