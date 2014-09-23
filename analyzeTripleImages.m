%%  load syncing data
dataFolder='E:\DATA\3DwormData\BrainScanner20140911_182843';
[bf2fluorIdx,fluorAll,bfAll]=YamlFlashAlign(dataFolder);
hiResData=highResTimeTraceAnalysis(dataFolder);


hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));
bfFlashTime=bfAll.frameTime(bfAll.flashLoc);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);
hi2lowTimeOffset=hiResFlashTime-bfFlashTime;
hiResData.frameTime=hiResData.frameTime-hi2lowTimeOffset(1);
hiResFlashTime=(hiResData.frameTime(hiResData.flashLoc));

fluor2bfTimeOffset=fluorFlashTime-bfFlashTime;
fluorAll.frameTime=fluorAll.frameTime-fluor2bfTimeOffset(1);
fluorFlashTime=fluorAll.frameTime(fluorAll.flashLoc);

%works amazingly well, flash frames off by less than .05 seconds. 

%% load centerline data
centerLineFile=dir([dataFolder filesep '*centerline*']);
centerLineFile={centerLineFile.name}';
if length(centerLineFile)>1
    centerlineFile=uipickfiles('FilterSpec',dataFolder);
    centerline=load(centerlineFile{1},'centerline');
    centerline=centerline.centerline;
else
centerline=load([dataFolder filesep centerLineFile],'centerline');
centerline=centerline.centerline;
end
%% load alignment data
display('Select Low Res Alignment')
lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load(lowResFluor2BF{1});
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);

Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowRes=load(Hi2LowRes{1});
t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');

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
    movies=uipickfiles;
    behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

behaviorVidObj = VideoReader(behaviorMovie);
 fluorVidObj= VideoReader(fluorMovie);

Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

%% make lookup tables for indices
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'PCHIP');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));


firstFullFrame=find(~isnan(bfIdxLookup),1,'first');
firstFullFrame=max(firstFullFrame,find(~isnan(fluorIdxLookup),1,'first'));

%%
stretchSize=15;
lowResFolder=[dataFolder filesep 'lowResFolder'];
hiResActivityFolder=[dataFolder filesep 'hiResActivityFolder'];
hiResSegmentFolder=[dataFolder filesep  'hiResSegmentFolder'];
p=.001;
W=ones(1,100);
W([1,2,3,98,99,100])=10;

mkdir(lowResFolder);
mkdir(hiResActivityFolder);
mkdir(hiResSegmentFolder);
frames=firstFullFrame:length(hiResData.frameTime); %start 1750, 12000 good too, 13000 for 3d
frames(ismember(frames,hiResData.flashLoc))=[];
movieFlag=0;
imH=NaN(1,4);
lineH=NaN(1,4);
if movieFlag
vidOut=VideoWriter([dataFolder filesep 'HiMagOnly.avi']);
vidOut.FrameRate=20;

open(vidOut);
end


for iFrame=frames
    %%
    iTime=hiResData.frameTime(iFrame);
    %interpolate using time to get low res idx 
    bfIdx=bfIdxLookup(iFrame);
fluorIdx=fluorIdxLookup(iFrame);
hiResIdx=hiResLookup(iFrame);
  status=fseek(Fid,2*hiResIdx*1024^2,-1);

  pixelValues=fread(Fid,1024^2,'uint16',0,'l');

hiResImage=(reshape(pixelValues,1024,1024)');
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
fluorFrame2=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
    'OutputView',lowResFluor2BF.Rsegment);
fluorFrame2=imwarp(fluorFrame2,Hi2LowRes.t_concord,  'OutputView',Hi2LowRes.Rsegment);
fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));

%  bfFrame=imwarp(fluorFrame,lowResFluor2BF.t_concord,...
%     'OutputView',lowResFluor2BF.Rsegment);
%%
CLcurrent=[interp2(squeeze(centerline(:,1,:))',1:100,repmat(bfIdx,1,100))',...
interp2(squeeze(centerline(:,2,:))',1:100,repmat(bfIdx,1,100))'];
CLcurrent=[csaps(1:100,CLcurrent(:,1),p,1:100,W)'...
    csaps(1:100,CLcurrent(:,2),p,1:100,W)'];
CLcurrent=[interp1(CLcurrent(:,1),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')',...
    interp1(CLcurrent(:,2),-stretchSize+1:100+stretchSize,'*PCHIP','extrap')'];

endpts=CLcurrent([1,length(CLcurrent)],:);
[maxPoint(1),maxPoint(2)]=find(fluorFrame==max(fluorFrame(:)),1,'first');
tip2maxDistance=pdist2(endpts,maxPoint);
if tip2maxDistance(2)<tip2maxDistance(1)
    CLcurrent=flipud(CLcurrent);
end


CLdistance=[0;cumsum(sqrt(sum(diff(CLcurrent).^2,2)))];
totLengthPix=max(CLdistance);
CLcurrent=[interp1(CLdistance,CLcurrent(:,1),1:1:totLengthPix)',...
    interp1(CLdistance,CLcurrent(:,2),1:1:totLengthPix)'];



[~,newX,newY]=wormStraightening(CLcurrent,[],40,1);
    [~,newXHi,newYHi]=wormStraightening(CLcurrent(1:min(300,length(CLcurrent)),:),[],40,.2);

% CLnew=([fliplr(CLcurrent) ones(length(CLcurrent),1)]*lowResBF2FluorT.T);
% CLnew=CLnew(:,[2,1]);

   [CLnewY,CLnewX]=transformPointsInverse(lowResFluor2BF.t_concord...
       , CLcurrent(:,2),CLcurrent(:,1));
   
   CLnew=[CLnewX CLnewY];
   [CLHighResY,CLHighResX]=transformPointsForward(Hi2LowRes.t_concord...
       , CLcurrent(:,2),CLcurrent(:,1));
   
   CLHighRes=[CLHighResX, CLHighResY];
   
   CLHighRes=CLHighRes(CLHighResX<1024 & CLHighResY<1024 & ...
       CLHighResX>0 & CLHighResY>0,:);
   
   [ newYHi,newXHi]=transformPointsForward(Hi2LowRes.t_concord...
       , newYHi,newXHi);
   newYHi=newYHi-rect1(1)-1;
   newXHi=newXHi-rect1(2)-1;
%    n
%       [ newYHiA,newXHiA]=transformPointsInverse(S2AHiRes.t_concord...
%        , newYHi,newXHi);
   %%
   highResInterp=interp2(segmentChannel,newYHi,newXHi);
   highResActivityInterp=interp2(activityChannel,newYHi,newXHi);
      lowResFluorHiInterp=interp2(double(fluorFrame2),newYHi,newXHi);
      lowResFluorHiInterp=pedistalSubtract(lowResFluorHiInterp,5);
      lowResFluorHiInterp(isnan(lowResFluorHiInterp))=0;


%[maxY,maxX]=find(lowResFluorHiInterp==max(lowResFluorHiInterp(:)),1,'first');
%chan vesse..for some reason
% margin = 10;
% smooth_weight = .0001; 
% image_weight = 1000; 
% delta_t = 4; 
% if ~exist('phi','var')
% phi = zeros(size(lowResFluorHiInterp)); 
% phi(margin:end-margin, margin:end-margin) = 1; 
% phi = ac_reinit(phi-.5); 
% 
%     phi = ac_ChanVese_model(max(phi(:))*normalizeRange(lowResFluorHiInterp),...
%         phi, smooth_weight, image_weight, delta_t, 20); 
% else
%       phi = ac_ChanVese_model(max(phi(:))*normalizeRange(lowResFluorHiInterp),...
%         phi, smooth_weight, image_weight, delta_t, 10);  
% end

%  [C,h]=contour(phi,[0,0],'w');axis equal;
  fluorThresh=lowResFluorHiInterp>(max(lowResFluorHiInterp(:))/2);
  if any(fluorThresh(:))
stats=regionprops(fluorThresh,lowResFluorHiInterp,'Area','WeightedCentroid');
centroid=stats([stats.Area]==max([stats.Area])).WeightedCentroid;
 maxX=round(centroid(1));
 maxY=round(centroid(2));
  end
  
% maxY=find(sum(lowResFluorHiInterp,2)==max(sum(lowResFluorHiInterp,2)));
%  maxX=find(sum(lowResFluorHiInterp,1)==max(sum(lowResFluorHiInterp,1)));

lowResFluorHiInterp=rectCrop(lowResFluorHiInterp,[maxX-150,maxY-200,maxX+150, maxY+200]);
highResInterp=rectCrop(highResInterp,[maxX-150,maxY-200,maxX+150, maxY+200]);
highResActivityInterp=rectCrop(highResActivityInterp,[maxX-150,maxY-200,maxX+150, maxY+200]);

jointIm=[1+max(highResActivityInterp(:))/60*normalizeRange(highResActivityInterp),...
    max(highResInterp(:))/1700*normalizeRange(highResInterp)];
jointIm(isnan(jointIm))=0;
jointIm(jointIm==1)=0;
hiResImage=max(highResInterp(:))/1700*normalizeRange(hiResImage);
%% plotting 4 images 
% subplot(2,2,1);
%     if (ishandle(imH(1)) && ishandle(lineH(1)));
%                 
% unfreezeColors
% set(imH(1),'Cdata',hiResImage);
% caxis([0,1])
%      %   set(lineH(1),'Xdata',CLHighRes(:,2),'Ydata',CLHighRes(:,1));
%         colormap hot
%         freezeColors
%     else
%     
%     imH(1)=imagesc(hiResImage);
%     caxis([0,1]);
%      %   axis equal;axis off;
% axis off; axis equal
%     hold on
%    % lineH(1)=plot(CLHighRes(:,2),CLHighRes(:,1),'g');
%         hold off
% colormap hot
%     freezeColors(subplot(2,2,1))
% 
%     end
%     title([ num2str( hiResData.frameTime(iFrame),'%.2f ') ' s']...
%         ,'fontsize',20,'color','w');
%     subplot(2,2,3)
%         if ishandle(imH(2));
% unfreezeColors
% set(imH(2), 'Cdata',jointIm);
% colormap hot
%    colormap(interp1([hot;circshift(hot,[0,1])],1:2:128));
%    caxis([0,2]);
% freezeColors
% 
%         else
%             
%     imH(2)=imagesc(jointIm);
%    colormap(interp1([hot;circshift(hot,[0,1])],1:2:128));
%       caxis([0,2]);
% freezeColors
%  %  freezeColors(subplot(2,2,3))
%     axis equal
%     axis off
%         end
% %         
%          title([ num2str( round(25*(hiResData.Z(iFrame)+1))) ' um  '...
%              num2str( hiResData.frameTime(iFrame),'%.2f ') ' s'] ...
%             ,'fontsize',20,'color','w')
%     subplot(2,2,2);
%         if (ishandle(imH(3)) && ishandle(lineH(3)));
% unfreezeColors
% set(imH(3),'Cdata',fluorFrame);
% set(lineH(3),'Xdata',CLnew(:,2),'Ydata',CLnew(:,1));
% colormap hot
% freezeColors
%         else
%     imH(3)=imagesc(fluorFrame);
%         axis equal;colormap hot; hold on;axis off
%     lineH(3)=plot(squeeze(CLnew(:,2)),squeeze(CLnew(:,1)),'g');
%     hold off
%     freezeColors(subplot(2,2,2));
%         end
%  %   title(num2str(fluorIdx));
%     subplot(2,2,4);
%     if  (ishandle(imH(4)) && ishandle(lineH(4)))
%         unfreezeColors
% set(imH(4),'Cdata',bfFrame);
% set(lineH(4),'Xdata',CLcurrent(:,2),'Ydata',CLcurrent(:,1));
% colormap gray
% freezeColors
%     else
%     imH(4)=imagesc(bfFrame);
%     
%         axis off;axis equal;   hold on
%     lineH(4)=plot(squeeze(CLcurrent(:,2)),squeeze(CLcurrent(:,1)),'g');
% %plot(newY',newX','black')
%     hold off
%     colormap gray
%     freezeColors(subplot(2,2,4))
%     spaceplots
 %    set(gcf,'color','black');
%         drawnow;
% 
%     end
    
    %%
    if movieFlag
    vidFrame=getframe(gcf);
    writeVideo(vidOut,vidFrame);
    end
    %    title(num2str(bfIdx));

    fileName=['image' num2str(iFrame,'%3.5d') '.tif'];
    tiffwrite([lowResFolder filesep fileName],lowResFluorHiInterp,'tif',0);
    tiffwrite([hiResActivityFolder filesep fileName],highResActivityInterp,'tif',0);
    tiffwrite([hiResSegmentFolder filesep fileName],highResInterp,'tif',0);
    
    
end
close(vidOut);