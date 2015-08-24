% worm=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResSegmentFolder3Dtest_raw\image00500.tif');
% activity=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResActivityFolder3Dtest_raw\image00500.tif');

% load folder with data
 dataFolder=uipickfiles;
 dataFolder=dataFolder{1};
%dataFolder='F:\20141212\BrainScanner20141212_145951\';
%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
%% load video and alignemtn data



zindexer=@(x,s) x./(s)+1;


rows=1200;
cols=600;
nPix=rows*cols;


[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);


    vidInfo.bfAll=bfAll;
    vidInfo.fluorAll=fluorAll;
    vidInfo.hiResData=hiResData;
    

%% make time alignments
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear','extrap');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear','extrap');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));



stack2BFidx=bfIdxLookup(diff(hiResData.stackIdx)==1);
stack2fluoridx=fluorIdxLookup(diff(hiResData.stackIdx)==1);

%stack2BFidx=stack2BFidx(~isnan(stack2BFidx))
BF2stackIdx=interp1(stack2BFidx,1:max(hiResData.stackIdx),bfIdxList,'nearest');
fluor2stackIdx=interp1(stack2fluoridx,1:max(hiResData.stackIdx),fluorIdxList,'nearest');


%% load centerline
behaviorFolder=dir([dataFolder filesep 'Behavior*']);
behaviorFolder=behaviorFolder([behaviorFolder.isdir]);
behaviorFolder=[dataFolder filesep behaviorFolder(1).name];
centerlineFile=dir([behaviorFolder filesep 'center*']);
centerlineFile=[behaviorFolder filesep centerlineFile.name];
centerline=load(centerlineFile);
CLfieldNames=fieldnames(centerline);
CLfieldIdx=cellfun(@(x) ~isempty(strfind(x,'centerline')),CLfieldNames);
CLoffsetIdx=cellfun(@(x) ~isempty(strfind(x,'off')),CLfieldNames);
if any(CLoffsetIdx)
CLoffset=centerline.(CLfieldNames{CLoffsetIdx});


else
   CLoffset=0; 
end

%%

minStart=max([min(BF2stackIdx(CLoffset+1)) min(fluor2stackIdx)])+1;

%% load alignment data
try
    alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
catch
display('Select Low Res Alignment')

lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
%lowResFluor2BF=load('Y:\CommunalCode\3dbrain\registration\20141212LowResBehavior2Fluor.mat');
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);


display('Select Hi to Low Fluor Res Alignment')
Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load(Hi2LowResF{1});
%Hi2LowResF=load('Y:\CommunalCode\3dbrain\registration\20141212HighResS2LowResFluorBeads.mat');


% display('Select Hi to Low Res Alignment')
% 
% Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
% Hi2LowRes=load(Hi2LowRes{1});
% t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');
 display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
%S2AHiRes=load('Y:\CommunalCode\3dbrain\registration\20141212HiResS2A.mat');
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;

%%


%%
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;

save([dataFolder filesep 'alignments'],'alignments');

end

%% 
fiducialFile=uipickfiles('filterspec',dataFolder);
fiducialPoints=load(fiducialFile{1});
fiducialPoints=fiducialPoints.fiducialPoints;


%% calculate offset between frame position and zposition
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%%

startStack=minStart;
endStack=max(hiResData.stackIdx);
destination= 'CLstraight_20150813';
imageFolder2=[dataFolder filesep destination];
mkdir(imageFolder2);
show=0;
stackRange= startStack:endStack;

pointStats=repmat(struct(),1,length(stackRange));

%% do first image in range
tic
show=1;

counter=50;
[V,pointStatsOut,Vtemplate,side,lastOffset,Vbw]=...
    WormCLStraighten_6(dataFolder,destination,vidInfo,...
    alignments,[],[],zOffset,minStart+counter,[],[],show);
poinStatsFields={'straightPoints','rawPoints','stackIdx','pointIdx',...
    'Rintensities','Volume','controlPoints'};
%poinStatsFields=fieldnames(pointStatsOut);
for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=pointStatsOut.(field);

    
    
end

show=0;




display(['Finished image ' num2str(startStack,'%3.5d') ' in ' num2str(toc) 's'])
[~,clusterFolder]=fileparts(dataFolder);
clusterFolder=['/scratch/tmp/jnguyen/' clusterFolder];
save([dataFolder filesep 'startWorkspace'])

%% par for through all other frames
%subfiducialPoints=fiducialPoints(stackRange);
%parforprogress(length(stackRange)-1);
if ~isfield(pointStats(1),'stackIdx')
    missingIdx=ones(1,length(pointStats));
else
missingIdx=cellfun(@(x) isempty(x),{pointStats.stackIdx})';
end
%stackRange2=stackRange(missingIdx);
parfor counter=1:500%length(stackRange);
    if  missingIdx(counter)
    
 %   parforprogress
%progressbar((iStack-startStack)/(endStack-startStack));
             iStack=stackRange(counter);
display(['Starting'  num2str(iStack,'%3.5d') ])
     try
          tic
%change indexing for better parfor 
         ctrlPoints=[];%subfiducialPoints{counter};
[~,pointStatsOut,~,~]=...
    WormCLStraighten_6(dataFolder,destination,vidInfo,...
    alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,show);

for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=pointStatsOut.(field);

    
end




display(['Finished image ' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])


     catch ME
        ME
        display(['Error in Frame' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])

    end
    else
       iStack=stackRange(counter);

    %   display(['Frame ' num2str(iStack,'%3.5d') ' already done!'])

    end
    
end
save([imageFolder2 filesep 'PointsStats'],'pointStats');

%%



