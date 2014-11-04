Ztype='traingle';
imSize=[1200,600];
row=imSize(1);
col=imSize(2);
%%  load syncing data
mostRecent=getappdata(0,'mostRecent');
dataFolder=uipickfiles('FilterSpec',mostRecent);
dataFolder=dataFolder{1};

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,imSize(1),imSize(2));
end
z2ImageIdxOffset=-11;


%% load alignment data
display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;



%% prep .dat file




%%
%progressbar(0,0)
stackIdx=hiResData.stackIdx;
mkdir([dataFolder filesep 'stackDataWhole']);

%%

options.thresh1=.03; %initial Threshold
options.hthresh=-.001; %threshold for trace of hessian.
options.minObjSize=140;
options.maxObjSize=Inf;
options.watershedFilter=1;
options.filterSize=[5 5 3];
options.pad=9;
options.noise=1;
options.show=0;
options.maxSplit=1;
options.minSphericity=.55;
options.valleyRatio=.8;
options.scaleFactor=[1,1,3];
spikeBuffer=0;
meanSliceFiltLevel=3;

gaussianFilter=fspecial('gaussian',[30,30],5);
gaussianFilter=convnfft(gaussianFilter,permute(gausswin(6,2),[2,3,1]));
gaussianFilter=gaussianFilter/sum(gaussianFilter(:));


%%
mkdir([dataFolder filesep 'stackDataWhole'])
outputData=[];
outputData.centroids=[];
outputData.Rintensities=[];
outputData.Gintensities=[];
outputData.Volume=[];
outputData.time=[];
outputData.wormMask=[];
outputDat.zPlaneIdx=[];
outputData.zMax=[];
outputData.zMin=[];
output0=outputData;
overwriteFlag=1;
stackList=2:max(stackIdx)-1;

if ~overwriteFlag;
    d=dir([dataFolder filesep 'stackDataWhole' filesep 'stack*']);
    d={d.name}';
    oldIdx=cell2mat(cellfun(@(x) str2double(x(6:9)),d,'UniformOutput',false));
    stackList=stackList((~ismember(stackList,oldIdx)'));
end
groupSize=50;
%%
for iGroup=2:length(stackList)/groupSize
    
    group=groupSize*(iGroup-1):groupSize*(iGroup);
    
    group=group(group<length(stackList));
    group=group(group>1);
    
    subStackList=stackList(group);
    outputData=repmat(output0,length(subStackList),1);
    
    parfor i=1:length(group);
        
        try
            
            iStack=subStackList(i);
            tic
            %%
            imageIdx=find(stackIdx==iStack);
            imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
            worm=[];activity=[];
            zVoltage=hiResData.Z(imageIdx);
            time=hiResData.frameTime(imageIdx);
            Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

            for iSlice=1:length(imageIdx)
                %read proper image file
                hiResIdx=imageIdx(iSlice)+z2ImageIdxOffset;
                status=fseek(Fid,2*hiResIdx*row*col,-1);
                pixelValues=fread(Fid,row*col,'uint16',0,'l');
                hiResImage=(reshape(pixelValues,row,col));
                %hiResImage=hiResImage-backgroundImage;
                %hiResImage(hiResImage<0)=0;
                
                %align 2 halves of dv2
                activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
                segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
                activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
                activityChannel=(activityChannel);
                %build stack
                worm(:,:,iSlice)=segmentChannel;
                activity(:,:,iSlice)=activityChannel;
            end
            fclose(Fid);
            worm=pedistalSubtract(worm,5);
            activity=pedistalSubtract(activity,5);
            
            
            %% do segmentation
            wormMask=WormSegmentHessian3d_rescale(worm,options);
            %%
            %look up intensities on both channels
            wormLabelMask=bwlabeln(wormMask,6);
            wormcc=bwconncomp(wormMask);
            stats=regionprops(wormcc,worm,'WeightedCentroid','Area','PixelIdxList');
            centroids=reshape([stats.WeightedCentroid],3,[])';
            
            %smooth out activity for interpolating, this is equivalent to expanding the
            %ROI's.
            activity=convnfft(activity,gaussianFilter,'same');
            
            
            Rintensities=cellfun(@(x) mean(worm(x)),[wormcc.PixelIdxList])';
            Gintensities=cellfun(@(x) mean(activity(x)),[wormcc.PixelIdxList])';
            
            %interpolate Z properly and scale
            realTime=interp1(time,centroids(:,3)); %arb scaling for now
            zPlaneIdx=centroids(:,3);
            %[~,~,zPix]=cellfun(@(x) ind2sub(size(K),x),{stats.PixelIdxList}','Uniformoutput',false);
            %zMax=cellfun(@(x) max(x), zPix);
            %zMin=cellfun(@(x) min(x), zPix);
            centroids(:,3)=50*(1+interp1(zVoltage,centroids(:,3))); %arb scaling for now
            Volume=[stats.Area]';
            
            %save outputs in unique file
            outputFile=[dataFolder filesep 'stackDataWhole' filesep 'stack' num2str(iStack,'%04d') 'data'];
            outputData(i).centroids=centroids;
            outputData(i).Rintensities=Rintensities;
            outputData(i).Gintensities=Gintensities;
            outputData(i).Volume=Volume;
            ouptutData(i).centroidTime=realTime;            
            outputData(i).wormMask=wormLabelMask;

            outputData(i).metaData.time=time;
            outputData(i).metaData.zPlaneIdx=zPlaneIdx;
            outputData(i).metaData.zVoltage=zVoltage;
            outputData(i).metaData.iFrame=imageIdx;
            %outputData(i).zMax=zMax;
            %outputData(i).zMin=zMin;
            %parsavestruct(outputFile,outputData(i));
            display(['Completed stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);
            
        catch ME
            display(['Error in stack' num2str(iStack,'%04d') 'in ' num2str(toc) ' seconds']);
            ME
        end
        
    end
    %%
    for i=1:length(outputData);
        if ~isempty(outputData(i).centroids)
            iStack=subStackList(i);
            outputFile=[dataFolder filesep 'stackDataWhole' filesep 'stack' num2str(iStack,'%04d') 'data'];
            parsavestruct(outputFile,outputData(i));
            display(['Saving' num2str(iStack,'%04d')]);
            
        end
    end
    
end



%%


