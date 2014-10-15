function runTrack(imFolder,metaDataFolder)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minDist=30;
params.good=20;
params.mem=10;
%imFolder=getappdata(0,'imFolder');
if nargin==1;
    metaDataFolder=[];
end


matFiles=dir([imFolder filesep '*.mat']);
%setappdata(0,'matFiles',matFiles);
%imFiles=dir([imFolder filesep '*.tif']);
%setappdata(0,'imFiles',imFiles);


%%
trackData=[];
trackIdx=0;
progressbar(0);

for imat=1:1:length(matFiles)
    trackIdx=trackIdx+1;
    matName=matFiles(imat).name;
    load([imFolder filesep matName]);
    if ~isempty(metaDataFolder);
        metaData=load([metaDataFolder filesep 'image0' matName(6:9)]);
        metaData=metaData.metaData;
        zPlaneIdx=zPlaneIdx-metaData.midPlane;
        zVoltage=metaData.zVoltage-interp1(metaData.zVoltage,metaData.midPlane);
        newZ=(1+interp1(zVoltage,zPlaneIdx))*50;
        centroids(:,3)=newZ;

    end
    
    

    tracks=[centroids,Gintensities,Rintensities,trackIdx*ones(size(Gintensities))];
    params.dim=size(centroids,2);

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
ndim=size(centroids,2);
%progressbar(0);
for iTrack=1:nTracks
 %   progressbar(iTrack/nTracks)
 
    t=trackOutput((trackOutput(:,end)==iTrack),end-1);
    centroid=trackOutput((trackOutput(:,end)==iTrack),1:ndim);
    green=trackOutput((trackOutput(:,end)==iTrack),ndim+1);
    red=trackOutput((trackOutput(:,end)==iTrack),ndim+2);
    
    cellOutput(iTrack).time=t;
    cellOutput(iTrack).centroid=centroid;
    cellOutput(iTrack).green=green;
    cellOutput(iTrack).red=red;
    subplot(2,1,1);
    plot(t,cellOutput(iTrack).green)
    subplot(2,1,2);
        plot(t,cellOutput(iTrack).red)
        
drawnow
 %   pause(1)

end
%%
save([imFolder filesep 'trackOutput'],'trackOutput','cellOutput')


