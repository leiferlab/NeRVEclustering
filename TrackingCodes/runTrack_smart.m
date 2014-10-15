function runTrack_smart(imFolder,metaDataFolder)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nargin==1;
    metaDataFolder=[];
end
%%
minDist=15;
params.good=20;
params.mem=10;
%imFolder=getappdata(0,'imFolder');


matFiles=dir([imFolder filesep 'stack*.mat']);
%setappdata(0,'matFiles',matFiles);
%imFiles=dir([imFolder filesep '*.tif']);
%setappdata(0,'imFiles',imFiles);


%%
trackData1=[];
trackData=[];
idxAll=[];
trackIdx=0;
progressbar(0);
endTime=155;%length(matFiles);

matData=load([imFolder filesep matFiles(1).name]);
centroids1=matData.centroids;
idx1=3*ones(size(centroids1,1),1);
idx1(centroids1(:,2)<300)=2;
idx1(centroids1(:,2)<200)=1;

 
 %%
 
for imat=2:1:endTime
    
    
    trackIdx=trackIdx+1;
    matName=matFiles(imat).name;
    load([imFolder filesep matName]);
%     if ~isempty(metaDataFolder);
%         metaData=load([metaDataFolder filesep 'image0' matName(6:9)]);
%         metaData=metaData.metaData;
%         if isfield(metaData,'metaData');
%             metaData=metaData.metaData;
%         end
%         
%         zPlaneIdx=zPlaneIdx-metaData.midPlane;
%         zVoltage=metaData.zVoltage-interp1(metaData.zVoltage,metaData.midPlane);
%         newZ=(1+interp1(zVoltage,zPlaneIdx))*50;
%         centroids(:,3)=newZ;
% 
%     end
%     
   

     idx=3*ones(size(centroids,1),1);
          idx(centroids(:,2)<300)=2;

     idx(centroids(:,2)<200)=1;
     
     offset=-mean(centroids(idx==2,1:2))+[150,250];
          centroids(idx==2,1:2)=bsxfun(@plus,centroids(idx==2,1:2),offset);

     %    c=lines(max(idx));
%    
%     c=c(idx,:);    
%        


    tracks=[centroids,Gintensities,Rintensities,(1:length(centroids))',trackIdx*ones(size(Gintensities))];
    
    params.dim=size(centroids,2);

    trackData=[trackData;tracks];
    idxAll=[idxAll;idx];
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
for iIdx=1:max(idxAll);
for i=1:10
    try
        trackOutput1{iIdx}=track(trackData(idxAll==iIdx,:),minDist*(1-.1*i),params);
       if iIdx>1;
        trackOutput1{iIdx}(:,end)= trackOutput1{iIdx}(:,end)+trackOutput1{iIdx-1}(end,end)+1;
    end    
        break
    catch
        display(['reducing minDist by factor of ' num2str(i)]);
    end


end
end
%%


trackOutput=cell2mat(trackOutput1');
% trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
% 
% badtracks=find(trackLengths<minTrack);
% badtracks=any(bsxfun(@eq, trackOutput(:,end),badtracks'),2);
% 
% trackOutput(badtracks,:)=[];
%  trackLengths=accumarray(trackOutput(:,end),ones(size(trackOutput(:,end))));
% [ trackIdx,ia,ib]=unique(trackOutput(:,end));
% trackOutput(:,end)=ib;
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
save([imFolder filesep 'trackOutput2'],'trackOutput','cellOutput')


