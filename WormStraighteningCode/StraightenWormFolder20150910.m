% worm=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResSegmentFolder3Dtest_raw\image00500.tif');
% activity=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResActivityFolder3Dtest_raw\image00500.tif');

% load folder with data
 dataFolder=uipickfiles;
 dataFolder=dataFolder{1};
%dataFolder='F:\20141212\BrainScanner20141212_145951\';%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
%% load video and alignemtn data

rows=1200;
cols=600;
nPix=rows*cols;


    [bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);

%bundle timing data

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

topLimit=min(max(hiResData.stackIdx),length(stack2BFidx));
BF2stackIdx=interp1(stack2BFidx,1:topLimit,bfIdxList,'nearest');
fluor2stackIdx=interp1(stack2fluoridx,1:topLimit,fluorIdxList,'nearest');


%% load centerline
[centerline, CLoffset]=loadCLBehavior(dataFolder);

%%

minStart=max([min(BF2stackIdx(CLoffset+1)) min(fluor2stackIdx)])+1;

%% load alignment data
try
    %try to load the alignment file in the folder, otherwise, select them
    %individual in the registration folder
    alignments=load([dataFolder filesep 'alignments']);
alignments=alignments.alignments;
catch
display('Select Low Res Alignment')

lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);

display('Select Hi to Low Fluor Res Alignment')
Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load(Hi2LowResF{1});

 display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;


%% if you select them individually, bundle them and save it for later use
alignments.lowResFluor2BF=lowResFluor2BF;
alignments.S2AHiRes=S2AHiRes;
alignments.Hi2LowResF=Hi2LowResF;

save([dataFolder filesep 'alignments'],'alignments');

end
%%
% If there is a hand annotated set for comparison, select here, otherwise
%skip
if 0
fiducialFile=uipickfiles('filterspec',dataFolder);
fiducialPoints=load(fiducialFile{1});
fiducialPoints=fiducialPoints.fiducialPoints;
end

%% calculate offset between frame position and zposition
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%% Create straightened destination folder
destination= ['CLstraight_' datestr(now,'yyyymmdd')];
imageFolder2=[dataFolder filesep destination];
mkdir(imageFolder2);

%% Select range of volumes to analyze

startStack=minStart;
endStack=max(hiResData.stackIdx);
stackRange= startStack:endStack;
pointStats=combinePointStatsFiles(imageFolder2,length(stackRange));
%% do a sample image in range, for observation and for reference
tic
show=1;

counter=1; %which volume to do

%Run straighten and segmentation on one volume
[V,pointStatsOut,Vtemplate,side,lastOffset,Vbw]=...
    WormCLStraighten_10(dataFolder,destination,vidInfo,...
    alignments,[],[],zOffset,minStart+counter,[],[],show);


% fields in the pointstats file that will be used
poinStatsFields={'straightPoints','rawPoints','stackIdx','pointIdx',...
    'Rintensities','Volume','controlPoints'};

for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=pointStatsOut.(field);
end

show=0;
display(['Finished image ' num2str(startStack,'%3.5d') ' in ' num2str(toc) 's'])
%save initial workspace for testing and later use
save([dataFolder filesep 'startWorkspace'])

%% possibly used for cluster straightening at a later point,
[~,clusterFolder]=fileparts(dataFolder);
clusterFolder=['/scratch/tmp/jnguyen/' clusterFolder];

%% par for to analyze all other frames
overWrite=0; % set overwrite flag
% find all the missing pointstats
if ~isfield(pointStats(1),'stackIdx')
    missingIdx=ones(1,length(pointStats));
else
missingIdx=cellfun(@(x) isempty(x),{pointStats.stackIdx})';
end


parfor counter=1:length(stackRange);
    tic
    %run on on frames that are not already done in the pointStats variable
    if  missingIdx(counter)
             iStack=stackRange(counter);
        psfileName=[dataFolder filesep destination filesep 'pointStats'...
            num2str(iStack,'%3.5d') '.mat'];
        %if ps file has already been saved, and we are not overwriting it,
        %load it. Otherwise, run wormstraightening
        if exist(psfileName,'file') && ~overWrite
            pointStatsOut=load(psfileName);
            pointStatsOut=pointStatsOut.pointStats;
           display(['Frame ' num2str(iStack,'%3.5d') ' already done!'])

        else

            display(['Starting'  num2str(iStack,'%3.5d') ])
            try
                ctrlPoints=[];%subfiducialPoints{counter};
                [~,pointStatsOut,~,~]=...
                    WormCLStraighten_10(dataFolder,destination,vidInfo,...
                    alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,show);
            catch ME
                ME
                display(['Error in Frame' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])
                pointStatsOut=pointStats(counter);
            end
            
        end
    
        if length(pointStatsOut.straightPoints)<250
for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=pointStatsOut.(field);
end
        end



display(['Finished image ' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])


    else
       iStack=stackRange(counter);

       display(['Frame ' num2str(iStack,'%3.5d') ' already done!'])

    end
    
end
save([dataFolder filesep 'PointsStats'],'pointStats');

%%



for counter=1:length(pointStats)
if length(pointStats(counter).straightPoints)>250
    for iFields=1:length(poinStatsFields)
    field=poinStatsFields{iFields};
    pointStats(counter).(field)=[];
    end
end
end
    