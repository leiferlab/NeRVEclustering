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

    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,rows,cols);




%% make time alignments

[hiImageIdx,ib]=unique(hiResData.imageIdx);

hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));






%%


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


%% calculate offset between frame position and zposition
zWave=hiResData.Z;
zWave=gradient(zWave);
zWave=smooth(zWave,10);
[ZSTDcorrplot,lags]=(crosscorr(abs(zWave),hiResData.imSTD,40));
ZSTDcorrplot=smooth(ZSTDcorrplot,3);
zOffset=lags(ZSTDcorrplot==max(ZSTDcorrplot));

%%

startStack=1;
endStack=max(hiResData.stackIdx);
destination= 'CLstraight_20160124fix';
imageFolder2=[dataFolder filesep destination];
mkdir(imageFolder2);
show=0;
stackRange= startStack:endStack;
pointStats=combinePointStatsFiles(imageFolder2,length(stackRange));

%% do first im+age in range
tic
show=1;

counter=150;
[V,pointStatsOut,Vtemplate,side,Vbw]=...
    WormCLStraighten_noCL(dataFolder,destination,hiResData,...
    alignments,[],[],zOffset,counter,[],[],show);
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
overWrite=1;
lastOffset=[0 0];

if ~isfield(pointStats(1),'stackIdx') || overWrite
    missingIdx=ones(1,length(pointStats));
else
missingIdx=cellfun(@(x) isempty(x),{pointStats.stackIdx})';
end
%stackRange2=stackRange(missingIdx);
parfor counter=1:length(stackRange);
    tic
    if  missingIdx(counter)
             iStack=stackRange(counter);
        psfileName=[dataFolder filesep destination filesep 'pointStats'...
            num2str(iStack,'%3.5d') '.mat'];
        if exist(psfileName,'file') && ~overWrite
            pointStatsOut=load(psfileName);
            pointStatsOut=pointStatsOut.pointStats;
           display(['Frame ' num2str(iStack,'%3.5d') ' already done!'])

        else
            %   parforprogressscan
            %progressbar((iStack-startStack)/(endStack-startStack));
            display(['Starting'  num2str(iStack,'%3.5d') ])
            try
                %change indexing for better parfor
                ctrlPoints=[];%subfiducialPoints{counter};
                [~,pointStatsOut]=...
                    WormCLStraighten_noCL_fx2(dataFolder,destination,hiResData,...
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
    