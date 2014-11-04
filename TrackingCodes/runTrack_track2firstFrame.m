function runTrack_track2firstFrame(imFolder,metaDataFolder)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nargin==1;
    metaDataFolder=[];
end
%%
minDist=15;

%imFolder=getappdata(0,'imFolder');


matFiles=dir([imFolder filesep '*.mat']);
%setappdata(0,'matFiles',matFiles);
%imFiles=dir([imFolder filesep '*.tif']);
%setappdata(0,'imFiles',imFiles);


%%
trackData1=[];
trackData=[];
idxAll=[];
trackIdx=0;
endTime=length(matFiles);

load([imFolder filesep matFiles(1).name]);
centroids1=centroids;
% idx1=3*ones(size(centroids1,1),1);
% idx1(centroids1(:,2)<300)=2;
% idx1(centroids1(:,2)<200)=1;

     tracks1=[centroids1,Gintensities,Rintensities,(1:length(centroids1))',ones(size(Gintensities))];
    params.dim=size(centroids1,2);
 params.quiet=1;
 %%
 trackOutput1=repmat({0},endTime,1);
 progressbar(0);

parfor imat=2:endTime
    try
    
    matName=matFiles(imat).name;
    stackData=load([imFolder filesep matName]);

centroids=stackData.centroids;
Rintensities=stackData.Rintensities;
Gintensities=stackData.Gintensities;
%[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(centroids, centroids1, 3, [1, 0.1, 0.01], [0.00008, 0.00008, 0.00008], [0 0 0], 10);
Transformed_M=centroids;
tracks=[Transformed_M,Gintensities,Rintensities,(1:length(centroids))',2*ones(size(Gintensities))];
    trackInput=cat(1,tracks1,tracks);
    for i=15:-1:1
    try
            trackOutput1{imat}=track(trackInput,i,params);
            break
    catch
    end
    end
    particleIdx=trackOutput1{imat}(:,end-1);
    particleIdx(particleIdx==2)=imat;
    trackOutput1{imat}(:,end-1)=particleIdx;
    display(['Finished stack: ' num2str(imat)]);
    catch
    end

end

%%

trackOutput2=trackOutput1;
firstTrack=trackOutput1{2};
trackTime=firstTrack(:,end-1);
particleIdx=firstTrack(:,end-2);
trackIdx=firstTrack(:,end);
particleIdx=particleIdx(trackTime==1);
[~,ib]=sort(particleIdx);
trackIdx=trackIdx(trackTime==1);
sortedTrackIdx=trackIdx(ib);
trackOutput=firstTrack;
%%
for imat=1:endTime
    try
    temp=trackOutput2{imat};
    tempTime=temp(:,end-1);
tempparticleIdx=temp(:,end-2);
tempTrackIdx=temp(:,end);
tempparticleIdx=tempparticleIdx(tempTime==1);
[~,ib]=sort(tempparticleIdx);
tempTrackIdx=tempTrackIdx(tempTime==1);
tempTrackIdx(ib)=sortedTrackIdx;
temp(tempTime==1,end)=tempTrackIdx;
time2=[false;diff(tempTime)>0];
time1=[diff(tempTime)>0;false];
temp(time2,end)=temp(time1,end);

trackOutput=cat(1,trackOutput,temp(time2,:));

if mod(imat,10)==0
    display(['Finished ' num2str(imat)]);
end

    catch
    end
    
end

trackOutput=sortrows(trackOutput,8);

%%
%save([imFolder filesep 'trackData'],'trackData');
%params.dim=size(centroids,2);
%%
centroids=[];

%%


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
%     subplot(2,1,1);
%     plot(t,cellOutput(iTrack).green)
%     subplot(2,1,2);
%         plot(t,cellOutput(iTrack).red)
%         
% drawnow
 %   pause(1)


end
%%
save([imFolder filesep 'trackOutput'],'trackOutput','cellOutput')


