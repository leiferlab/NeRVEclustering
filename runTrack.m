function runTrack(imFolder)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minDist=30;
params.good=100;
params.mem=30;
params.dim=2;
%imFolder=getappdata(0,'imFolder');

matFiles=dir([imFolder filesep '*.mat']);
%setappdata(0,'matFiles',matFiles);
%imFiles=dir([imFolder filesep '*.tif']);
%setappdata(0,'imFiles',imFiles);

trackData=[];
trackIdx=0;
progressbar(0);

for imat=1:1:length(matFiles)
    trackIdx=trackIdx+1;
    load([imFolder filesep matFiles(imat).name]);
    tracks=[centroids,Gintensities,Rintensities,trackIdx*ones(size(Gintensities))];
    trackData=[trackData;tracks];
    progressbar((imat)/length(matFiles));
    
%     figure(1)
%     if any(centroids)
%         imagesc(wormMask);
%         hold on
%     scatter(centroids(:,1),centroids(:,2),'g');
% 
%     drawnow
%     end
end
save([imFolder filesep 'trackData'],'trackData');
%params.dim=size(centroids,2);
%%
for i=1:10
    try
        trackOutput=track(trackData,minDist*(1-.1*i),params);
        break
    catch
        display(['reducing minDist by factor of ' num2str(i)]);
    end
end

% trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
% 
% badtracks=find(trackLengths<minTrack);
% badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
% 
% trackOutput(badtracks,:)=[];
%  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
[ trackIdx,ia,ib]=unique(trackOutput(:,end));
trackOutput(:,end)=ib;
%%
nTracks=max(trackOutput(:,end));
nTime=max(trackOutput(:,end-1));
%progressbar(0);
for iTrack=1:nTracks
 %   progressbar(iTrack/nTracks)
    t=trackOutput((trackOutput(:,end)==iTrack),end-1);
    centroid=trackOutput((trackOutput(:,end)==iTrack),1:2);
    green=trackOutput((trackOutput(:,end)==iTrack),3);
    red=trackOutput((trackOutput(:,end)==iTrack),4);
    
    cellOutput(iTrack).time=t;
    cellOutput(iTrack).centroid=centroid;
    cellOutput(iTrack).green=green;
    cellOutput(iTrack).red=red;
    subplot(2,1,1);
    plot(t,cellOutput(iTrack).green)
    subplot(2,1,2);
        plot(t,cellOutput(iTrack).red)
        

    cell
    pause(1)

end
%%
save([imFolder filesep 'trackOutput'],'trackOutput','cellOutput')


