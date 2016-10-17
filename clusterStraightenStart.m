
function clusterStraightenStart(dataFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates Initial files for Straightening images
% Inputs:
%   dataFolder- file path to data folder containing raw .dat file
%               timing results, behavior folder, lowmag folder,
%               alignments file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%which volume to do for initial, not too close to begining because start
%might not have all video feeds tracking properly
counter=300; 
%% bundle timing data
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder);
vidInfo.bfAll=bfAll;
vidInfo.fluorAll=fluorAll;
vidInfo.hiResData=hiResData;

%% make time alignments
bf_ft=bfAll.frameTime;
fluor_ft=fluorAll.frameTime;
hi_ft=hiResData.frameTime;

bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluor_ft);
bfIdxLookup=interp1(bf_ft,bfIdxList,hi_ft,'linear','extrap');
fluorIdxLookup=interp1(fluor_ft,fluorIdxList,hi_ft,'linear','extrap');
stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
stack2fluoridx=fluorIdxLookup(diff(hiResData.stackIdx)==1);

topLimit=min(max(hiResData.stackIdx),length(stack2BFidx));
BF2stackIdx=interp1(stack2BFidx,1:topLimit,bfIdxList,'nearest');
fluor2stackIdx=interp1(stack2fluoridx,1:topLimit,fluorIdxList,'nearest');


%% load find the first frame with every video in it
[~, CLoffset]=loadCLBehavior(dataFolder);
first_frame=max([min(BF2stackIdx(CLoffset+1)) min(fluor2stackIdx)])+1;
test_frame=first_frame+counter;

%% load alignment data

%try to load the alignment file in the folder, otherwise, select them
%individual in the registration folder
alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;

% deal with missing backgrounds
if ~isfield(alignments,'background');
    alignments.background=0;
    save([dataFolder filesep 'alignments'],'alignments');
end

%% calculate offset between frame position and zposition
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,100);
image_std=hiResData.imSTD;
image_std=image_std-mean(image_std);
image_std(image_std>100)=0;
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),image_std,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%% Create straightened destination folder
destination= ['CLstraight_' datestr(now,'yyyymmdd')];
imageFolder2=[dataFolder filesep destination];
mkdir(imageFolder2);

%% do a sample image in range, for observation and for reference

tic
%Run straighten and segmentation on one volume
[~,~,Vtemplate,side,~,~]=WormCLStraighten_11(...
    dataFolder,destination,vidInfo,alignments,[],zOffset,test_frame,'left',0); 
        
%we need these of the outputs to save for the cluster straightening

display(['Finished image ' num2str(test_frame,'%3.5d') ' in ' num2str(toc) 's'])
%save initial workspace from first sample for later use
save([dataFolder filesep 'startWorkspace'],...
    'destination', 'Vtemplate', 'zOffset', 'side','vidInfo')


