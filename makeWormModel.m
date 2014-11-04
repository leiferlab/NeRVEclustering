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
z2ImageIdxOffset=-7;


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
%%
contourFile=dir([dataFolder filesep '*ontours*']);
contourFile={contourFile.name}';
if length(contourFile)~=1
    contourFile=uipickfiles('FilterSpec',dataFolder);
    contours=load(contourFile{1});
    contours=contours.contours;
else
    contours=load([dataFolder filesep contourFile{1}]);
    contours=contours.contours;
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
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,bfAll.frameTime,'linear');

firstFullFrame=find(~isnan(fluorIdxLookup),1,'first');
firstFullFrame=firstFullFrame+1100;
%%
stretchSize=25;
%lowResFolder=[dataFolder filesep 'lowResFolder'];
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

%%
testSize=0;
counter=1;
Cout=[];
imagesAll=[];
for iFrame=firstFullFrame:testSize+firstFullFrame;
    tic

        %%
    try 
        bfIdx=iFrame+30*counter;
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
        
        %%
   CLcurrent=centerline(:,:,bfIdx);
  %          CLcurrent=squeeze(trimmean(CLcurrent,20,1));
                      CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
            CLdistance=1+CLdistance/CLdistance(end)*99;
            CLcurrent=[interp1(CLdistance,CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')',...
                interp1(CLdistance,CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')'];
      
            %%
            %  progressbar(iStack/max(stackIdx),iSlice/length(imageIdx))
              fluorIdx=round(fluorIdxLookup(bfIdx));
            Ctemp=contours(fluorIdx).C;
            Ctemp(Ctemp(:,1)==0,:)=[];
            %  hiResImage2=imwarp(hiResImage,Hi2LowRes.t_concord,...
            %     'OutputView',Hi2LowRes.Rsegment);
            
            fluorFrame=read(fluorVidObj,round(fluorIdx));
            bfFrame = read(behaviorVidObj,round(bfIdx));
            fluorFrame=fluorFrame(:,:,1);
            bfFrame=bfFrame(:,:,1);
            fluorFrame2=imwarp(fluorFrame,Hi2LowResF.t_concord,...
                'OutputView',Hi2LowResF.Rsegment);
            fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
             
%                bfFrame=imwarp(bfFrame,Hi2LowRes.t_concord,...
%                    'OutputView',Hi2LowRes.Rsegment);
%                bfFrame=bfFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
%%

            % CLcurrent=[csaps(1:100,CLcurrent(:,1),p,1:100,W)'...
            %     csaps(1:100,CLcurrent(:,2),p,1:100,W)'];
  %    CLcurrent=smooth2a(CLcurrent,10,0);
            endpts=CLcurrent([1,length(CLcurrent)],:);
            [maxPointx,maxPointy]=find(fluorFrame==max(fluorFrame(:)),1,'first');
            maxPoint=[maxPointx,maxPointy];
            tip2maxDistance=pdist2(endpts,maxPoint);
            CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
            totLengthPix=max(CLdistance);
            CLcurrent=[interp1(CLdistance,CLcurrent(:,1),1:1:totLengthPix)',...
                interp1(CLdistance,CLcurrent(:,2),1:1:totLengthPix)'];
            
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
            %optional correction to hi res coordinates
            CLHighRes=bsxfun(@plus,CLHighRes,hiResCorrection);
               [ newYHiF,newXHiF]=transformPointsForward(Hi2LowResF.t_concord...
                ,newYHi,newXHi);
            
            [ newYHiFlo,newXHiFlo]=transformPointsInverse(Hi2LowResF.t_concord...
                , newYHiF,newXHiF);  
            [CtempY,CtempX]=transformPointsForward(Hi2LowResF.t_concord,...
                Ctemp(:,1),Ctemp(:,2));
            
            %cut out coordinates into rect for segment
            CtempY=CtempY-rect1(1)+hiResCorrection(2)-1;
            CtempX=CtempX-rect1(2)+hiResCorrection(1)-1;

            newYHiF=newYHiF-rect1(1)+hiResCorrection(2)-1;
            newXHiF=newXHiF-rect1(2)+hiResCorrection(1)-1;
    
            %       [ newYHiA,newXHiA]=transformPointsInverse(S2AHiRes.t_concord...
            %        , newYHi,newXHi);
            %%
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
                maxX=400;%round(centroid(1));
               % end
                maxY=round(centroid(2));
            end
            % maxY=find(sum(lowResFluorHiInterp,2)==max(sum(lowResFluorHiInterp,2)));
            %  maxX=find(sum(lowResFluorHiInterp,1)==max(sum(lowResFluorHiInterp,1)));
            
            lowResFluorHiInterp=rectCrop(lowResFluorHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);
            bfHiInterp=rectCrop(bfHiInterp,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow]);

            %%
            newYHiF=rectCrop(newYHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
            newXHiF=rectCrop(newXHiF,[maxX-xWindow,maxY-yWindow,maxX+xWindow, maxY+yWindow],nan);
            k=dsearchn([newXHiF(:),newYHiF(:)],[CtempX,CtempY]);
[CtempX,CtempY]=ind2sub(size(newXHiF),k);

Ctemp=[CtempX,CtempY];

%[Transformed_M, multilevel_ctrl_pts, multilevel_param] = gmmreg_L2_multilevel(model, C2, 2, [.1, 0.01], [0.000008, 0.00008, 0.0000008], [0 0],1,0);

%compile data for each image
        imagesAll{counter}=lowResFluorHiInterp;

            Cout{counter}=Ctemp;
            counter=counter+1;
     
        imagesc(lowResFluorHiInterp)
        
        hold on
        scatter(CtempY,CtempX,'xr')
        drawnow
hold off
            
        display(['Finished frame: ' num2str(iFrame) ' in ' num2str(toc) 's']);
    catch me
        display(['Error frame: ' num2str(iFrame)]);
    end
    
end

refImage=lowResFluorHiInterp;
save([dataFolder filesep 'lowResFluorContour'],'Cout','refImage')
