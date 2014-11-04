%%analyzes triple images with just the centerline from the low mag data






Ztype='traingle';
imSize=[1200,600];
row=imSize(1);
col=imSize(2);
nPix=row*col;
%%  load syncing data
dataFolder='E:\DATA\3DwormData\20141020\20141020_1610Data';

[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);

if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle3(dataFolder,imSize(1),imSize(2));
end

hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
[~,bf2Hi]=flashTimeAlign2(bfFlashTime,hiResFlashTime);
f_hiResTime=fit(hiResFlashTime,bfFlashTime(bf2Hi),'poly1');
hiResData.frameTime=f_hiResTime(hiResData.frameTime);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

[~,bf2fluor]=flashTimeAlign2(bfFlashTime,fluorFlashTime);
f_fluorTime=fit(fluorFlashTime,bfFlashTime(bf2fluor),'poly1');
fluorAll.frameTime=f_fluorTime(fluorAll.frameTime);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

%works amazingly well, flash frames off by less than .05 seconds.
% NVM, only worked for early data set, has problem


%% load centerline data
centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile={centerLineFile.name}';
if length(centerLineFile)>1
    centerlineFile=uipickfiles('FilterSpec',dataFolder);
    centerline=load(centerlineFile{1},'centerline');
    centerline=centerline.centerline;
else
    centerline=load([dataFolder filesep centerLineFile{1}],'centerline');
    centerline=centerline.centerline;
end
%% load alignment data
display('Select Low Res Alignment')
lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);
display('Select Hi to Low Res Alignment')

Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowRes=load(Hi2LowRes{1});
t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');
display('Select Hi Res Alignment')

S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load(S2AHiRes{1});
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;

%% load background image
try
    backgroundImageFile=dir([dataFolder filesep '*backgroundImage*']);
    if isempty(backgroundImageFile) || length(backgroundImageFile)>1
        display('Select Background Mat File')
        backgroundImageFile=uipickfiles('FilterSpec',dataFolder);
        backgroundImageFile=backgroundImageFile{1};
    else
        backgroundImageFile=[dataFolder filesep backgroundImageFile.name];
    end
    backgroundImage=load(backgroundImageFile);
    backgroundImage=backgroundImage.backgroundImage;
catch
    backgroundImage=0;
end




%% prep vidobj for avi files and .dat file

%search in dataFolder for avi's without the HUDS, one with the fluor and
%one without.
aviFiles=dir([dataFolder filesep '*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
if length(aviFiles)==2
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    behaviorMovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluorMovie=[dataFolder filesep aviFiles{aviFluorIdx}];
else
    display('Select avi files, behavior and then low mag fluor');
    movies=uipickfiles('FilterSpec',dataFolder);
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);

Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

%% make lookup tables for indices
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));

%%
stretchSize=25;
%lowResFolder=[dataFolder filesep 'lowResFolder'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder3D'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder3D'];
lookupFolderX=[dataFolder filesep  'lookupFolderX'];
lookupFolderY=[dataFolder filesep  'lookupFolderY'];
Ralign=0;
metaFolder=[dataFolder filesep  'metaDataFolder3D'];
mkdir(metaFolder);
mkdir(hiResActivityFolder);
mkdir(hiResSegmentFolder);
mkdir(lookupFolderX);
mkdir(lookupFolderY);
frames=1750:length(hiResData.frameTime); %start 1750, 12000 good too, 13000 for 3d
frames(ismember(frames,hiResData.flashLoc))=[];
movieFlag=0;
plotFlag=0;
saveFlag=1;
lookupSave=1;
imH=NaN(1,4);
lineH=NaN(1,4);
if movieFlag
    plotFlag=1;
    vidOut=VideoWriter([dataFolder filesep 'HiMagOnly.avi']);
    vidOut.FrameRate=20;
    
    open(vidOut);
end
spikeBuffer=0;
meanSliceFiltLevel=4;
yWindow=300;
xWindow=200;
hiResCorrection=[0,0];
splitFlag=0;
stackIdx=hiResData.stackIdx;

gaussianKernal=fspecial('gaussian', 80,20);

%%
%progressbar(0,0)
metaData=[];
metaData.metaData=[];
stackIdx2=interp1(stackIdx,hiResLookup);
Zlookup=interp1(hiResData.Z,hiResLookup);
timeLookup=interp1(hiResData.frameTime,hiResLookup);
metaData=repmat(metaData,1,max(stackIdx2));
parfor iStack=2:max(stackIdx);
    try
        tic
        imageIdx=find(stackIdx==iStack);
        imageIdx=imageIdx(ismember(imageIdx,hiResData.imageIdx));
        imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
        worm=[];
        zPos=Zlookup(imageIdx);
        time=timeLookup(imageIdx);
        zPos=hiResData.Z(imageIdx);
        time=hiResData.frameTime(imageIdx);
        
        activity=[];
        highResInterpStack2=[];
        lowResFluorHiInterpStack=[];
        highResInterpStack=[];
        hiResMean=[];
        newYHiFout=[];
        highResActivityInterpStack=[];
        newXHiFout=[];
        segmentStack=[];
        activityStack=[];
        Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);
        
        %%
        for iSlice=1:length(imageIdx)
            %  progressbar(iStack/max(stackIdx),iSlice/length(imageIdx))
            iFrame=imageIdx(iSlice);
            iTime=hiResData.frameTime(iFrame);
            %interpolate using time to get low res idx
            hiResIdx=(iFrame);
            
            bfIdx=round(bfIdxLookup(hiResIdx));
            fluorIdx=round(fluorIdxLookup(hiResIdx));
            status=fseek(Fid,2*(hiResIdx)*nPix,-1);
            
            pixelValues=fread(Fid,nPix,'uint16',0,'l');
            
            hiResImage=(reshape(pixelValues,row,col));
            hiResImage=hiResImage-backgroundImage;
            hiResImage(hiResImage<0)=0;
            activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
            segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
            activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
            activityChannel=pedistalSubtract(activityChannel,5);
            segmentStack(:,:,iSlice)=segmentChannel;
            activityStack(:,:,iSlice)=activityChannel;
            %
            %  hiResImage2=imwarp(hiResImage,Hi2LowRes.t_concord,...
            %     'OutputView',Hi2LowRes.Rsegment);
            
            fluorFrame=read(fluorVidObj,round(fluorIdx));
            bfFrame = read(behaviorVidObj,round(bfIdx));
            fluorFrame=fluorFrame(:,:,1);
            bfFrame=bfFrame(:,:,1);
            
            
            fluorFrame2=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
                'OutputView',lowResFluor2BF.Rsegment);
            fluorFrame2=imwarp(fluorFrame2,Hi2LowRes.t_concord,  'OutputView',Hi2LowRes.Rsegment);
            fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
            
            %  bfFrame=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
            %     'OutputView',lowResFluor2BF.Rsegment);
            %%
            CLcurrent=[interp2(squeeze(centerline(:,1,:))',1:100,repmat(bfIdx,1,100))',...
                interp2(squeeze(centerline(:,2,:))',1:100,repmat(bfIdx,1,100))'];
            % CLcurrent=[csaps(1:100,CLcurrent(:,1),p,1:100,W)'...
            %     csaps(1:100,CLcurrent(:,2),p,1:100,W)'];
            CLcurrent=[interp1(CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')',...
                interp1(CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')'];
       %     CLcurrent=smooth2a(CLcurrent,10,0);
            endpts=CLcurrent([1,length(CLcurrent)],:);
            [maxPointx,maxPointy]=find(fluorFrame==max(fluorFrame(:)),1,'first');
            maxPoint=[maxPointx,maxPointy];
            tip2maxDistance=pdist2(endpts,maxPoint);
            if tip2maxDistance(2)<tip2maxDistance(1)
                CLcurrent=flipud(CLcurrent);
            end
            
            
            CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
            totLengthPix=max(CLdistance);
            CLcurrent=[interp1(CLdistance,CLcurrent(:,1),1:1:totLengthPix)',...
                interp1(CLdistance,CLcurrent(:,2),1:1:totLengthPix)'];
            
            
            
            [~,newX,newY]=wormStraightening(CLcurrent,[],80,10);
            
            [~,newXHi,newYHi]=wormStraightening(CLcurrent(1:min(Inf,length(CLcurrent)),:),[],80,10);
            [Xgrid, Ygrid]=meshgrid(1:.02:size(newXHi,2),1:.02:size(newXHi,1));
            newXHi=interp2(newXHi,Xgrid,Ygrid);
            newYHi=interp2(newYHi,Xgrid,Ygrid);
            
            % CLnew=([fliplr(CLcurrent) ones(length(CLcurrent),1)]*lowResBF2FluorT.T);
            % CLnew=CLnew(:,[2,1]);
            
            [CLnewY,CLnewX]=transformPointsInverse(lowResFluor2BF.t_concord...
                , CLcurrent(:,2),CLcurrent(:,1));
            
            CLnew=[CLnewX CLnewY];
            [CLHighResY,CLHighResX]=transformPointsForward(Hi2LowRes.t_concord...
                , CLcurrent(:,2),CLcurrent(:,1));
            
            CLHighRes=[CLHighResX, CLHighResY];
            CLHighRes=bsxfun(@plus,CLHighRes,hiResCorrection);
            
            CLHighRes=CLHighRes(CLHighResX<1024 & CLHighResY<1024 & ...
                CLHighResX>0 & CLHighResY>512,:);
            
            [ newYHiF,newXHiF]=transformPointsForward(Hi2LowRes.t_concord...
                ,newYHi,newXHi);
            
            newYHiF=newYHiF-rect1(1)+hiResCorrection(2)-1;
            newXHiF=newXHiF-rect1(2)+hiResCorrection(1)-1;
            %    n
            %       [ newYHiA,newXHiA]=transformPointsInverse(S2AHiRes.t_concord...
            %        , newYHi,newXHi);
            %%
            highResInterp=interp2(segmentChannel,newYHiF,newXHiF);
            highResActivityInterp=interp2(activityChannel,newYHiF,newXHiF);
            
            [ newYHiFlo,newXHiFlo]=transformPointsInverse(lowResFluor2BF.t_concord...
                , newYHi,newXHi);
            
            lowResFluorHiInterp=interp2(double(fluorFrame),newYHiFlo,newXHiFlo);
            lowResFluorHiInterp=pedistalSubtract(lowResFluorHiInterp,5);
            lowResFluorHiInterp(isnan(lowResFluorHiInterp))=0;
            
            
            
            fluorThresh=lowResFluorHiInterp>(max(lowResFluorHiInterp(:))/2);
            if any(fluorThresh(:))
                stats=regionprops(fluorThresh,lowResFluorHiInterp,'Area','WeightedCentroid');
                centroid=stats([stats.Area]==max([stats.Area])).WeightedCentroid;
                maxX=round(centroid(1));
                maxY=round(centroid(2));
            end
            % maxY=find(sum(lowResFluorHiInterp,2)==max(sum(lowResFluorHiInterp,2)));
            %  maxX=find(sum(lowResFluorHiInterp,1)==max(sum(lowResFluorHiInterp,1)));
            
            lowResFluorHiInterp=rectCrop(lowResFluorHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            highResInterp=rectCrop(highResInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            highResActivityInterp=rectCrop(highResActivityInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            newYHiFout(:,:,iSlice)=rectCrop(newYHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
            newXHiFout(:,:,iSlice)=rectCrop(newXHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
            
            
            highResInterp=pedistalSubtract(highResInterp);
            lowResFluorHiInterpStack(:,:,iSlice)=pedistalSubtract(lowResFluorHiInterp);
            highResInterpStack(:,:,iSlice)=highResInterp;
            highResActivityInterpStack(:,:,iSlice)=pedistalSubtract(highResActivityInterp);
            hiResMean(iSlice)=nanmean((highResInterp(:)));
            
        end
        fclose(Fid);
        
        
        
        
        %%
        stackSize=size(highResInterpStack,3);
        
        %find middle, ie where the ventral nerve cord is, its in the lower half of
        %the worm and is a peak in the I(z) plot, but not too close to top and
        %bottom
        
        if Ralign==1;
            tempStack=highResInterpStack;
        else
            tempStack=highResActivityInterpStack;
        end
        
        maxSearch=squeeze(nansum(nansum(tempStack,1),2));
        [~,iPeak]=max(maxSearch);
        midSearch=smooth(squeeze(nansum(nansum(tempStack(400:end,:,:),1),2)),4);
        [pVals,pLocs]=findpeaks(midSearch);
        pVals(pLocs<(iPeak-4) |pLocs>(iPeak+4))=0;
        [~,pmax]=max(pVals);
        if ~isempty(pmax)
            midPlane=pLocs(pmax(1));
        else
            midPlane=iPeak;
        end
        
        
        zProjWormBot=nanmax(tempStack(:,:,1:midPlane),[],3);
        zProjWormTop=nanmax(tempStack(:,:,midPlane:end),[],3);
        
        zProjWormMid=nanmax(tempStack(:,:,midPlane-3:midPlane+3),[],3);
        zProjWormMid(isnan(zProjWormMid))=0;
        zProjWormTop(isnan(zProjWormTop))=0;
        zProjWormBot(isnan(zProjWormBot))=0;


%         zProjWormMid=normalizeRange(smooth2a(zProjWormMid,10,30));
%         zProjWormTop=normalizeRange(smooth2a(zProjWormTop,10,40));
%         zProjWormBot=normalizeRange(smooth2a(zProjWormBot,10,40));
         zProjWormMid=normalizeRange(imfilter(zProjWormMid,gaussianKernal));
         zProjWormTop=normalizeRange(imfilter(zProjWormTop,gaussianKernal));
         zProjWormBot=normalizeRange(imfilter(zProjWormBot,gaussianKernal));
        
        
        %posterior is same for top and bottom
        PostProj=normalizeRange(zProjWormMid(300:end,:,:));
        zProjWormTop(300:end,:,:)=PostProj;
        zProjWormBot(300:end,:,:)=PostProj;
        
        Pbot2=ActiveContourFit2_1D(zProjWormTop,[],[]);
        Ptop2=ActiveContourFit2_1D(zProjWormBot,[],Pbot2);
        
        Iprofile=(smooth(interp2(smooth2a(nanmean(tempStack,3),50,50),Ptop2(:,1),Ptop2(:,2)),10));
        [~,Icenter]=max(Iprofile);
        
        Pbot2=interp1(Pbot2, Icenter-299:Icenter+300,'*linear','extrap');
        Ptop2=interp1(Ptop2, Icenter-299:Icenter+300,'*linear','extrap');
        
        %try to fix bias of contour to ventral nerve cord
        xVec=1:600;
        xCorrection=100./(1+exp(-(xVec'-333)/50));
        Pbot2(:,1)=Pbot2(:,1)-xCorrection;
        Ptop2(:,1)=Ptop2(:,1)-xCorrection;
        
%         Pbot2(301:end,1)=Pbot2(301:end,1)-(30-30*exp(-xVec/100))';
%         Ptop2(301:end,1)=Ptop2(301:end,1)-(30-30*exp(-xVec/100))';
%         
%         
        
        [~,newXbot,newYbot]=wormStraightening((fliplr(Pbot2)),[],150,1);
        [~,newXtop,newYtop]=wormStraightening((fliplr(Ptop2)),[],150,1);
        
        %%
        midPlaneBuffer=4;
        w=((1:stackSize)-midPlane+midPlaneBuffer)/(midPlaneBuffer*2);
        w=max(w,0);
        w=min(w,1);
        w=reshape(w,1,1,[]);
        [~,~,stackGrid]=meshgrid(1:size(newXtop,2),1:size(newXtop,1),1:stackSize);
        
        newXi=bsxfun(@times,newXbot,(1-w))+bsxfun(@times,newXtop,(w));
        newYi=bsxfun(@times,newYbot,(1-w))+bsxfun(@times,newYtop,(w));
        newXi2=interp3(newXHiFout,newYi,newXi,stackGrid);
        newYi2=interp3(newYHiFout,newYi,newXi,stackGrid);
        segmentStack2=interp3(segmentStack,newYi2,newXi2,stackGrid);
        activityStack2=interp3(activityStack,newYi2,newXi2,stackGrid);
        newYi2=uint16(newYi2);
        newXi2=uint16(newXi2);
        
        %   Iprofile=(smooth(max(smooth2a(sum(segmentStack2,3),10,10),[],2),30));
        %   [ymax,imax,ymin,imin]=extrema(Iprofile);
        %   imin=sort(imin);=
        %
        %   dips=[find(imin<imax(1),1,'last'),find(imin>imax(1),1,'first')];
        %   dips=imin(dips);
        %
        % p=[1:100,linspace(101,dips(1),100),linspace(dips(1)+1,dips(2),200),dips(2)+1:dips(2)+200];
        % segmentStack3=bsxfun(@(A,B) interp1(A,B),segmentStack2,p');
        
        
        %%
        if movieFlag
            vidFrame=getframe(gcf);
            writeVideo(vidOut,vidFrame);
        end
        %%
        
        %    title(num2str(bfIdx));
        metaData(iStack).metaData.time=time;
        metaData(iStack).metaData.zVoltage=zPos;
        metaData(iStack).metaData.iFrame=imageIdx;
        metaData(iStack).metaData.midPlane=midPlane;
        fileName=['image' num2str(iStack,'%3.5d') '.tif'];
        if saveFlag
            
            tiffwrite([hiResActivityFolder filesep fileName],single(activityStack2),'tif',0);
            tiffwrite([hiResSegmentFolder filesep fileName],single(segmentStack2),'tif',0);
            if lookupSave
                tiffwrite([lookupFolderX filesep fileName],newXi2,'tif',0);
                tiffwrite([lookupFolderY filesep fileName],newYi2,'tif',0);
                
            end
            
        end
        display(['Finished frame: ' num2str(iStack) ' in ' num2str(toc) 's']);
    catch me
        display(['Error frame: ' num2str(iStack)]);
    end
    
end

for iStack=1:max(stackIdx)
    matName=['image' num2str(iStack,'%3.5d')];
    parsave([metaFolder filesep matName],metaData(iStack),'metaData');
end


close(vidOut);