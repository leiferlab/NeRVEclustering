%%analyzes triple images with just the centerline from the low mag data



Ztype='traingle';
imSize=[1200,600];
row=imSize(1);
col=imSize(2);
nPix=row*col;
%%  load syncing data
dataFolder=uipickfiles;
dataFolder=dataFolder{1};

%% algin all videos,    CHECK Z2IMAGEIDXOFFSET
[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
z2ImageIdxOffset=-9;


%% load centerline data
centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile={centerLineFile.name}';
if length(centerLineFile)>1
            display('Select centerLine file');

    centerlineFile=uipickfiles('FilterSpec',dataFolder);
    centerline=load(centerlineFile{1},'centerline');
    centerline=centerline.centerline;
else
    centerline=load([dataFolder filesep centerLineFile{1}],'centerline');
    centerline=centerline.centerline;
end
%%
contourFile=dir([dataFolder filesep '*ontours*']);
contourFile={contourFile.name}';
if length(contourFile)~=1
    display('Select contour file');
    contourFile=uipickfiles('FilterSpec',dataFolder);
    contours=load(contourFile{1});
%    contours=contours.contours;
else
    contours=load([dataFolder filesep contourFile{1}]);
  %  contours=contours.contours;
end
%%
modelFile=dir([dataFolder filesep '*model*']);
modelFile={modelFile.name}';
if length(modelFile)~=1
        display('Select model file');

    modelFile=uipickfiles('FilterSpec',dataFolder);
    cModel=load(contourFile{1});
else
    cModel=load([dataFolder filesep modelFile{1}]);
    cModel=cModel.cModel;
end


%% load alignment data
display('Select Low Res Alignment')
lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);


display('Select Hi to Low Fluor Res Alignment')
Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load(Hi2LowResF{1});


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
%%
load([dataFolder filesep 'RefStack.mat'])

refstack=normalizeRange(refstack);
maxRef=max(pedistalSubtract(refstack),[],3);
maxRef=normalizeRange(maxRef);
maxRef(isnan(maxRef))=0;
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
hiResActivityFolder=[dataFolder filesep 'test' filesep 'hiResActivityFolder3Dtest'];
hiResSegmentFolder=[dataFolder filesep 'test' filesep  'hiResSegmentFolder3Dtest'];
lookupFolderX=[dataFolder filesep 'test' filesep  'lookupFolderXtest'];
lookupFolderY=[dataFolder filesep 'test' filesep  'lookupFolderYtest'];
metaFolder=[dataFolder filesep 'test' filesep  'metaDataFolder3Dtest'];
hiResSegmentFolderRaw=[dataFolder filesep 'test' filesep  'hiResSegmentFolder3Dtest_raw'];
hiResActivityFolderRaw=[dataFolder filesep 'test' filesep  'hiResActivityFolder3Dtest_raw'];
hiResSegmentFolderS1=[dataFolder filesep 'test' filesep  'hiResSegmentFolder3Dtest_S1'];
hiResActivityFolderS1=[dataFolder filesep 'test' filesep  'hiResActivityFolder3Dtest_S1'];
lowResFolder=[dataFolder filesep 'test' filesep  'lowResFluor'];
lowResBFFolder=[dataFolder filesep 'test' filesep  'lowResBF'];
mkdir(metaFolder);
mkdir(hiResActivityFolder);
mkdir(hiResSegmentFolder);
mkdir(lookupFolderX);
mkdir(lookupFolderY);
mkdir(hiResSegmentFolderRaw);
mkdir(hiResSegmentFolderS1);
mkdir(hiResActivityFolderRaw);
mkdir(hiResActivityFolderS1);
mkdir(lowResFolder);
mkdir(lowResBFFolder);

frames=1750:length(hiResData.frameTime); %start 1750, 12000 good too, 13000 for 3d
frames(ismember(frames,hiResData.flashLoc))=[];
movieFlag=0;
plotFlag=0;
saveFlag=1;
lookupSave=1;
imH=NaN(1,4);
lineH=NaN(1,4);
Ralign=2;
if movieFlag
    plotFlag=1;
    vidOut=VideoWriter([dataFolder filesep 'HiMagOnly.avi']);
    vidOut.FrameRate=20;
    
    open(vidOut);
end
spikeBuffer=0;
meanSliceFiltLevel=4;
yWindow=300;
xWindow=300;
hiResCorrection=[0,0];
splitFlag=0;
stackIdx=hiResData.stackIdx;

gaussianKernal=fspecial('gaussian', 80,20);
gaussianFilter=fspecial('gaussian',[100,100],30);
gaussianFilter=convnfft(gaussianFilter,permute(gausswin(10,3),[2,3,1]));

stackIdx2=interp1(stackIdx,hiResLookup);
Zlookup=interp1(hiResData.Z,hiResLookup);
timeLookup=interp1(hiResData.frameTime,hiResLookup);

%%  
Options.Registration='NonRigid';
  Options.Similarity='p';
  Options.MaxRef=20;
Options.Verbose=0;
Options.SigmaFluid=12;

%%
%progressbar(0,0)
interpMethod.method='nearest';
interpMethod.radius=40;
groupSize=80;% best to use factor of 8
for iGroup=1:max(stackIdx)/groupSize
metaData=[];
metaData.metaData=[];
group=groupSize*(iGroup-1):groupSize*(iGroup);

group=group(group<max(stackIdx));
group=group(group>1);

metaData=repmat(metaData,1,length(group));

parfor i=1:length(group);
    iStack=group(i);
    
    try
        %%
        tic
        imageIdx=find(stackIdx==iStack);
        imageIdx=imageIdx(ismember(imageIdx,hiResData.imageIdx));
        imageIdx=imageIdx(spikeBuffer+1:end-spikeBuffer);
        worm=[];
        zPos=Zlookup(imageIdx);
        time=timeLookup(imageIdx);
        zPos=hiResData.Z(imageIdx);
        time=hiResData.frameTime(imageIdx);
        maxXAll=[];maxYAll=[];
        activity=[];
        rawStack=[];
        highResInterpStack2=[];
        lowResFluorHiInterpStack=[];
        highResInterpStack=[];
        hiResMean=[];
        newYHiFout=[];
        lowResBFStack=[];
        highResActivityInterpStack=[];
        newXHiFout=[];
        segmentStack=[];
        activityStack=[];
        Cout=[];
        Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);
        
        %%
        CLIdx=round(bfIdxLookup(imageIdx));
        [xCLsearch,yCLsearch]=meshgrid(1:100,CLIdx);
               CLcurrent=cat(3,interp2(squeeze(centerline(:,1,:))',xCLsearch,yCLsearch),...
                interp2(squeeze(centerline(:,2,:))',xCLsearch,yCLsearch));
            CLcurrent=squeeze(trimmean(CLcurrent,20,1));
           CLcurrent=[interp1(CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')',...
                interp1(CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')'];
                  CLcurrent=[smooth(CLcurrent(:,1),5), smooth(CLcurrent(:,2),5)];
                        CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
          %  CLdistance=1+CLdistance/CLdistance(end)*99;
            CLcurrent=[interp1(CLdistance,CLcurrent(:,1),1:max(CLdistance),'*PCHIP','extrap')',...
                interp1(CLdistance,CLcurrent(:,2),1:max(CLdistance),'*PCHIP','extrap')'];
      
            %%
        for iSlice=1:length(imageIdx)
            %  progressbar(iStack/max(stackIdx),iSlice/length(imageIdx))
            iFrame=imageIdx(iSlice);
            iTime=hiResData.frameTime(iFrame);
            %interpolate using time to get low res idx
            hiResIdx=(iFrame);
            
            bfIdx=round(bfIdxLookup(hiResIdx));
            fluorIdx=round(fluorIdxLookup(hiResIdx));
         
            status=fseek(Fid,2*(hiResIdx+z2ImageIdxOffset)*nPix,-1);
            
            pixelValues=fread(Fid,nPix,'uint16',0,'l');
            
            hiResImage=(reshape(pixelValues,row,col));
            hiResImage=hiResImage-backgroundImage;
            hiResImage(hiResImage<0)=0;
            activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
            segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
            activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
            activityChannel=pedistalSubtract(activityChannel,5);
            %
            %  hiResImage2=imwarp(hiResImage,Hi2LowRes.t_concord,...
            %     'OutputView',Hi2LowRes.Rsegment);
            
            fluorFrame=read(fluorVidObj,round(fluorIdx));
            bfFrame = read(behaviorVidObj,round(bfIdx));
            fluorFrame=fluorFrame(:,:,1);
            bfFrame=bfFrame(:,:,1);
            fluorFrame2=imwarp(fluorFrame,Hi2LowResF.t_concord,...
                'OutputView',Hi2LowResF.Rsegment);
            fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
             
                   fluorbfFrame=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
                    'OutputView',lowResFluor2BF.Rsegment);
%                bfFrame=imwarp(bfFrame,Hi2LowRes.t_concord,...
%                    'OutputView',Hi2LowRes.Rsegment);
%                bfFrame=bfFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
%%


            %CLcurrent in BF reference
            
            [~,newX,newY]=wormStraightening(CLcurrent,[],80,10);
            [~,newXHi,newYHi]=wormStraightening(CLcurrent(1:min(Inf,length(CLcurrent)),:),[],80,10);
            [Xgrid, Ygrid]=meshgrid(1:.02:size(newXHi,2),1:.02:size(newXHi,1));
            newXHi=interp2(newXHi,Xgrid,Ygrid);
            newYHi=interp2(newYHi,Xgrid,Ygrid);
            
            %transform CL coordinates into fluorFrame
            [CLnewY,CLnewX]=transformPointsInverse(lowResFluor2BF.t_concord...
                , CLcurrent(:,2),CLcurrent(:,1));
            CLnew=[CLnewX CLnewY];
           [newYHi,newXHi]=transformPointsInverse(lowResFluor2BF.t_concord...
                , newYHi,newXHi);           
            
            
            %transform CL coordinates into HiRes Image
            [CLHighResY,CLHighResX]=transformPointsForward(Hi2LowResF.t_concord...
                ,CLnewY,CLnewX);
            CLHighRes=[CLHighResX, CLHighResY];
            %optional correction to hi res coordinates and cropping
            CLHighRes=bsxfun(@plus,CLHighRes,hiResCorrection-rect1([2,1])-1);

            
            %% transform cooridnate system
            [ newYHiF,newXHiF]=transformPointsForward(Hi2LowResF.t_concord...
                ,newYHi,newXHi);
            [ newYHiFlo,newXHiFlo]=transformPointsInverse(Hi2LowResF.t_concord...
                , newYHiF,newXHiF);  
            
            newYHiF=newYHiF-rect1(1)+hiResCorrection(2)-1;
            newXHiF=newXHiF-rect1(2)+hiResCorrection(1)-1;
    
            %       [ newYHiA,newXHiA]=transformPointsInverse(S2AHiRes.t_concord...
            %        , newYHi,newXHi);
            %%
            highResInterp=interp2(segmentChannel,newYHiF,newXHiF);
            highResActivityInterp=interp2(activityChannel,newYHiF,newXHiF);

            bfHiInterp=interp2(double(bfFrame),newYHi,newXHi);
            bfHiInterp=pedistalSubtract(bfHiInterp,5);
            bfHiInterp(isnan(bfHiInterp))=0;
     
            
            lowResFluorHiInterp=interp2(double(fluorFrame),newYHiFlo,newXHiFlo);
            lowResFluorHiInterp=pedistalSubtract(lowResFluorHiInterp,5);
            lowResFluorHiInterp(isnan(lowResFluorHiInterp))=0;
            
            fluorThresh=lowResFluorHiInterp>(max(lowResFluorHiInterp(:))/2);
            if any(fluorThresh(:)) 
                stats=regionprops(fluorThresh,lowResFluorHiInterp,'Area','WeightedCentroid');
                centroid=stats([stats.Area]==max([stats.Area])).WeightedCentroid;
              %  if iSlice==1
               %   maxX=round(centroid(1));
               % end
               maxX=400;
                maxY=round(centroid(2));
                maxXAll(iSlice)=maxX;
                maxYAll(iSlice)=maxY;
            end
            % maxY=find(sum(lowResFluorHiInterp,2)==max(sum(lowResFluorHiInterp,2)));
            %  maxX=find(sum(lowResFluorHiInterp,1)==max(sum(lowResFluorHiInterp,1)));
            
            lowResFluorHiInterp=rectCrop(lowResFluorHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            highResInterp=rectCrop(highResInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            highResActivityInterp=rectCrop(highResActivityInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            bfHiInterp=rectCrop(bfHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);

            %%
            newYHiF=rectCrop(newYHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
            newXHiF=rectCrop(newXHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);

% 
%   [Ireg,Bx,By,Fx,Fy] = register_images(lowResFluorHiInterp,refImage,Options);
%   Ireg2=movepixels(highResInterp,Fx,Fy);

% %%
% subSamp=randsample(length(Ctemp),length(Ctemp)/3);
% Transformed_M=Ctemp;
% %Transformed_M=bsxfun(@minus, Transformed_M,mean(Ctemp)-mean(cModel));
% [Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel( Transformed_M(subSamp,:),cModel, 3, [2, 0.01,.001], [0.0008, 0.00008, 0.0000008], [0 0 0],1,1);

%%
%[imgw, imgwr, map] = tpswarp(lowResFluorHiInterp, size(lowResFluorHiInterp), Ctemp(subSamp,:), Transformed_M, interpMethod)
 
 
%compile data for each image
            newYHiFout(:,:,iSlice)=newYHiF;
            newXHiFout(:,:,iSlice)=newXHiF;
            lowResBFStack(:,:,iSlice)=fluorFrame2;
            highResInterp=pedistalSubtract(highResInterp);
            lowResFluorHiInterpStack(:,:,iSlice)=pedistalSubtract(lowResFluorHiInterp);
            highResInterpStack(:,:,iSlice)=highResInterp;
            highResActivityInterpStack(:,:,iSlice)=pedistalSubtract(highResActivityInterp);
            hiResMean(iSlice)=nanmean((highResInterp(:)));
            
                        activityStack(:,:,iSlice)=activityChannel;
                        segmentStack(:,:,iSlice)=segmentChannel;
        end
        
        fclose(Fid);
        
      
 %       smoothedSegmentStack=smooth3(highResInterpStack,'gaussian',[101,101,11],50);
        %%
        stackSize=size(highResInterpStack,3);
        
        %find middle, ie where the ventral nerve cord is, its in the lower half of
        %the worm and is a peak in the I(z) plot, but not too close to top and
        %bottom
        
        if Ralign==1;
            tempStack=highResInterpStack;

        elseif Ralign==2;
            tempStack=sum(cat(4,highResActivityInterpStack,highResInterpStack),4);

        elseif Ralign==3
            tempStack=highResActivityInterpStack;
        else
            tempStack=lowResFluorHiInterpStack;
        end
  %%      
        for iIt=1:1;
        nanMap=isnan(tempStack);

        tempStack(isnan(tempStack))=0;         
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
        %%
        
        
        smoothedStack=convnfft(tempStack,gaussianFilter,'same');
        smoothedStack(nanMap)=nan;
        smoothedStack=colNanFill(smoothedStack);
smoothedCorr=nanmean(nanmean(smoothedStack,2),3);
smoothedCorr=max(smoothedCorr,1);
smoothedCorr=smooth(smoothedCorr,30);
smoothedCorr2=1-mean(mean(nanMap,2),3);
smoothedStack=normalizeRange(smoothedStack)+normalizeRange(bsxfun(@rdivide,smoothedStack,smoothedCorr));
    smoothedStack=normalizeRange(bsxfun(@times,smoothedStack,smoothedCorr2));

Pbot=ActiveContourFit2_psudo2D(normalizeRange(smoothedStack),[],[]);

        
        Iprofile=(smooth(interp2(smooth2a(nanmean(tempStack,3),50,50),Pbot(:,1,round(stackSize/2)),Pbot(:,2,round(stackSize/2))),3));
        [~,Icenter]=max(Iprofile);
        
        Pbot=interp1(Pbot, Icenter-299:Icenter+300,'*linear','extrap');
        
        %try to fix bias of contour to ventral nerve cord
        xVec=1:600;
       xCorrection=50./(1+exp(-(xVec'-333)/50));
       
       if lean(1)=='r'
       Pbot(:,1,:)=bsxfun(@minus,Pbot(:,1,:),xCorrection);
       elseif lean(1)=='l'
           Pbot(:,1,:)=bsxfun(@plus,Pbot(:,1,:),xCorrection);
       end
       
     
     
     
      %  Ptop2(:,1)=Ptop2(:,1)-xCorrection;
        
%         Pbot2(301:end,1)=Pbot2(301:end,1)-(30-30*exp(-xVec/100))';
%         Ptop2(301:end,1)=Ptop2(301:end,1)-(30-30*exp(-xVec/100))';
%         
%         
        

        %%
        [~,newXi,newYi]=wormStraightening(Pbot(:,[2,1],:),[],150,1);

        [~,~,stackGrid]=meshgrid(1:size(newXi,2),1:size(newXi,1),1:stackSize);
        
        [~,~,lineGrid]=meshgrid(1,1:size(Pbot,1),1:stackSize);
         reCL=[interp3(newXHiFout,Pbot(:,2,:),Pbot(:,1,:),lineGrid),...
             interp3(newYHiFout,Pbot(:,2,:),Pbot(:,1,:),lineGrid)];
 s=cat(1,zeros(1,1,stackSize),cumsum(sqrt(diff(reCL(:,1,:),[],1).^2+diff(reCL(:,2,:),[],1).^2),1));
 
 
        newXi2=interp3(newXHiFout,newYi,newXi,stackGrid);
        newYi2=interp3(newYHiFout,newYi,newXi,stackGrid);
        segmentStack2=interp3(segmentStack,newYi2,newXi2,stackGrid);
        activityStack2=interp3(activityStack,newYi2,newXi2,stackGrid);
        newYHiFout=(newYi2);
        newXHiFout=(newXi2);
        tempStack=pedistalSubtract(segmentStack2+activityStack2);
        end
        %%
        newYHiFout=uint16(newYHiFout);
        newXHiFout=uint16(newXHiFout);
        
        tempAlign=max(pedistalSubtract(segmentStack2),[],3);
        tempAlign=normalizeRange(tempAlign);
tempAlign(isnan(tempAlign))=0;
        
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
        metaData(i).metaData.time=time;
        metaData(i).metaData.zVoltage=zPos;
        metaData(i).metaData.iFrame=imageIdx;
        metaData(i).metaData.midPlane=midPlane;
        fileName=['image' num2str(iStack,'%3.5d') '.tif'];
        if saveFlag
            
            tiffwrite([hiResActivityFolder filesep fileName],single(activityStack2),'tif',0);
            tiffwrite([hiResSegmentFolder filesep fileName],single(segmentStack2),'tif',0);
            tiffwrite([hiResSegmentFolderRaw filesep fileName],single(segmentStack),'tif',0);
            tiffwrite([hiResSegmentFolderS1 filesep fileName],single(highResInterpStack),'tif',0);
            tiffwrite([hiResActivityFolderRaw filesep fileName],single(activityStack),'tif',0);
            tiffwrite([hiResActivityFolderS1 filesep fileName],single(highResActivityInterpStack),'tif',0);
            tiffwrite([lowResFolder filesep fileName],single(lowResFluorHiInterpStack),'tif',0);
           % tiffwrite([lowResBFFolder filesep fileName],single(lowResBFStack),'tif',0);
    
           
            if lookupSave
                
                
                [Ireg,Bx,By,Fx,Fy] = register_images(tempAlign,maxRef,Options);

                tiffwrite([lookupFolderX filesep fileName],single(Fx),'tif',0);
                tiffwrite([lookupFolderY filesep fileName],single(Fy),'tif',0);
                
            end
            
        end
        display(['Finished frame: ' num2str(iStack) ' in ' num2str(toc) 's']);
    catch me
        display(['Error frame: ' num2str(iStack)]);
    end
    
end

for iStack=1:length(metaData);
    
    matName=['image' num2str(group(iStack),'%3.5d')];
    parsave([metaFolder filesep matName],metaData(iStack),'metaData');
    display(['Saving ' matName]);

end

end
