function runTrack_tps_planeselect(imFolder,metaDataFolder, zRange)
% hObject    handle to runTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nargin==1;
    metaDataFolder=[];
end
%%
minDist=30;
zMin=zRange(1);
zMax=zRange(2);

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
endTime=length(matFiles);

load([imFolder filesep matFiles(2).name]);
centroids1=centroids;
% idx1=3*ones(size(centroids1,1),1);
% idx1(centroids1(:,2)<300)=2;
% idx1(centroids1(:,2)<200)=1;
centroids1(:,3)=zPlaneIdx;
zSelect=(zPlaneIdx>=zRange(1))&(zPlaneIdx<=zRange(2));
%centroids1=centroids1(zSelect,:);
tracks1=[centroids1,Gintensities,Rintensities,(1:size(centroids1,1))',ones(size(Gintensities))];
tracks1=tracks1(zSelect,:);
params.dim=size(centroids1,2);
params.good=2;
params.mem=2;
params1=params;
params1.mem=10;
params.quiet=1;
params.difficult=3e4;
params.excessive=4;
params1.good=2;
params1.excessive=4;



%% track each data set to first one (or 10 in this case)
trackOutput1=repmat({0},endTime,1);
trackData=tracks1;
initialNum=5;
if initialNum>1
for imat=2:initialNum;
    try
    matName=matFiles(imat).name;
    load([imFolder filesep matName]);

            if ~isempty(metaDataFolder);
            metaData=load([metaDataFolder filesep 'image0' matName(6:9)]);
            metaData=metaData.metaData;
            metaData=metaData.metaData;
            midPlane=metaData.midPlane;
            zVoltage=metaData.zVoltage;
            end
             coffset=median(centroids(centroids(:,2)>400,3));

            if isnan(coffset);
              coffset=median(centroids(:,3));
            end
            
         centroids(:,3)=zPlaneIdx;
zSelect=(zPlaneIdx>=zRange(1))&(zPlaneIdx<=zRange(2));
centroids=centroids(zSelect,:);
  %       centroids(:,3)=centroids(:,3)-50*(1+zVoltage(midPlane));
               %   zPlaneIdx=zPlaneIdx-metaData.midPlane;
         %  [ metaData.zVoltage(metaData.midPlane),interp1(metaData.zVoltage,metaData.midPlane)]

         %    zVoltage=metaData.zVoltage-interp1(metaData.zVoltage,metaData.midPlane);
           %  newZ=(1+interp1(zVoltage,zPlaneIdx))*50;
          %  centroids(:,3)=newZ;
            
       %     centroids(:,3)=centroids(:,3)-mean(centroids(:,3))+100;
Transformed_M=centroids;
[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(centroids, centroids1, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0 0], 1,0);
    tracks=[Transformed_M,Gintensities(zSelect),Rintensities(zSelect),find(zSelect),imat*ones(size(centroids(:,1)))];
    params.dim=size(centroids,2);
    trackData=[trackData;tracks];
    
    catch me
        display('Error caught')
        rethrow(me)
    end
% scatter3(centroids(:,1),centroids(:,2),centroids(:,3))
% drawnow
% hold on 
end
%%
trackOutput=nan;
counter=minDist;
while isnan(trackOutput)
    
        trackOutput=trackJN(trackData,counter,params1);
counter=counter-1;
        display(['reducing minDist, ' num2str(counter)]);
end
%%
[~,ib]=unique(trackOutput(:,end));
initialC=trackOutput(ib,:);
initialC(:,end-1)=1;
initialC2=initialC;
initialC=initialC(:,1:end-1);
initialC(:,end-1)=initialC2(:,end);
params.excessive=4;
centroids1=initialC(:,1:3);
else 
    initialC=tracks1;
end


%%
parfor imat=1:endTime
    try
    tic
    %%
    matName=matFiles(imat).name;
    stackData=load([imFolder filesep matName]);
    if isfield(stackData,'centroids')
        centroids=stackData.centroids;
        Rintensities=stackData.Rintensities;
        Gintensities=stackData.Gintensities;
        zPlaneIdx=stackData.zPlaneIdx;
    else
        centroids=[];
    end
    if ~isempty(centroids)
        
        if ~isempty(metaDataFolder);
%             metaData=load([metaDataFolder filesep 'image0' matName(6:9)]);
%             metaData=metaData.metaData;
%             zPlaneIdx=zPlaneIdx-metaData.midPlane;
%              zVoltage=metaData.zVoltage-interp1(metaData.zVoltage,metaData.midPlane);
%              newZ=(1+interp1(zVoltage,zPlaneIdx))*50;
%             centroids(:,3)=newZ;
        end
            coffset=median(centroids(centroids(:,2)>400,3));
            if isnan(coffset);
              coffset=nanmedian(centroids(centroids(:,2)>300,3));
            end
            
         centroids(:,3)=zPlaneIdx;
zSelect=(zPlaneIdx>=zRange(1))&(zPlaneIdx<=zRange(2));
centroids=centroids(zSelect,:);
        
%         
      %  centroids(:,3)=centroids(:,3)-mean(centroids(:,3))+100;

[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(centroids, centroids1, 3, [4, 0.2, 0.01], [0.0000008, 0.0000008, 0.0000008], [0 0 0], 1 ,0);
        T1=toc;
        tracks=[Transformed_M,Gintensities(zSelect),Rintensities(zSelect),find(zSelect),2*ones(nnz(zSelect),1)]
        trackInput=cat(1,initialC,tracks);
        tempOut=NaN;
        counter=20;
        while isnan(tempOut)    
            tempOut=trackJN(trackInput,counter,params);
            if isnan(tempOut)
            counter=counter-1;
            display('lowering dist to counter');

            if counter==1;
                break
            end
            end
            
        end
        trackOutput1{imat}=tempOut;
        trackOutput1{imat}(trackOutput1{imat}(:,end-1)==1,end-1)=0;
        trackOutput1{imat}(trackOutput1{imat}(:,end-1)==2,end-1)=imat;
        display(['Finished stack: ' num2str(imat) ' in ' num2str(toc) 's' ...
            ' Mindist: ' num2str(counter) 'Alignment done in :' num2str(T1) 's']);
    else
        trackOutput1{imat}=[];
        display(['Error stack: ' num2str(imat) ' in ' num2str(toc) 's']);
        
    end
    catch me
        me
                trackOutput1{imat}=[];
        display(['Error stack: ' num2str(imat) ' in ' num2str(toc) 's']);
        
    end
    
    
end

%% combine all tracks,

trackOutput2=trackOutput1;
firstTrack=initialC2; %take second track
trackTime=firstTrack(:,end-1);%get times
particleIdx=firstTrack(:,end-2); %get idx within
trackIdx=firstTrack(:,end);% get tracked idx
%particleIdx=particleIdx(trackTime==1);
[~,ib]=sort(particleIdx);
%trackIdx=trackIdx(trackTime==1);
sortedTrackIdx=trackIdx(ib); %match tracked idx with idx within
%% replace the index from track with initial identities
trackOutput=[];
for imat=1:endTime
    temp=trackOutput2{imat};
    if any(temp(:))
    tempTime=temp(:,end-1);
    tempparticleIdx=temp(:,end-2);
    tempTrackIdx=temp(:,end);
    tempparticleIdx=tempparticleIdx(tempTime==0);
    [tempparticleIdx,ib]=sort(tempparticleIdx);
    
    tempTrackIdx=tempTrackIdx(tempTime==0);
    tempTrackIdx(ib)=sortedTrackIdx(tempparticleIdx);
    temp(tempTime==0,end)=tempTrackIdx;
    time2=[false;diff(tempTime)>0];
    time1=[diff(tempTime)>0;false];
    temp(time2,end)=temp(time1,end);
    
    trackOutput=cat(1,trackOutput,temp(time2,:));
    
    if mod(imat,10)==0
        display(['Finished ' num2str(imat)]);
    end
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
ndim=size(centroids1,2);
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
       pause(2)
    
    
end
%%
save([imFolder filesep 'singlePlaneTrackOutput'],'trackOutput','cellOutput')


