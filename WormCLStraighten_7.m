
function [V,pointStats,Vproj,side,xyOffset2,wormBW2]=WormCLStraighten_7(dataFolder,destination,vidInfo,alignments,ctrlPoints,Vtemplate,zOffset,iStack,side,lastOffset,show)

%use different alignment than initial version, no need to crop after
%transformation
%added somethign to fix above

imageFolder2=[dataFolder filesep destination];
if isempty(lastOffset);
lastOffset=[0 0];
end
% if ~isempty(firstFrameParam)
% Vtemplate=firstFrameParam.Vproj
% 
% end
% 
% firstFrameParam.Vproj=Vproj;
% firstFrameParam.vRegion=vRegion;
% firstFrameParam.side=side;
% firstFrameParam.lastFrame=xyOffset2;

try
%% initial parameters
outputRadius=63.5;
outputRadiusZ=63.5;
outputLength=300;
CLsearchWindow=25;
zRatio=1/3;
zindexer=@(x,s) x./(s)+1;
    options.method='invdist';
    options.radius=20;
    options.power=1;
    options.thresh1=.05;
    options.minObjSize=50;
    options.maxObjSize=400;
    options.minSphericity=.80;
options.filterSize=[10 10 4];
options.power=1;
    options.prefilter=1;
    options.hthresh=0;
    
%% set up different kernals

 gaussKernal2=gausswin(200);
 gaussKernal2=convnfft(gaussKernal2,gaussKernal2');



Sfilter=max(gaussKernal2(:))-gaussKernal2;

Sfilter(Sfilter<.1)=-(.1-Sfilter(Sfilter<.1))*80;%Sfilter(Sfilter<.01)-.3;

Sfilter(Sfilter>.8)=0;
Sfilter(Sfilter>0)=1;%nnz(Sfilter<0)/nnz(Sfilter>0);
%Sfilter=Sfilter-mean(Sfilter(:));

Sfilter2=max(gaussKernal2(:))-gaussKernal2;
Sfilter2(Sfilter2<.01)=Sfilter2(Sfilter2<.01)-.3;
Sfilter2(Sfilter2>.6)=0;

Sfilter2(Sfilter2>0)=nnz(Sfilter2<0)/nnz(Sfilter2>0);
Sfilter2=Sfilter2-mean(Sfilter2(:));





%% recover alignments
lowResFluor2BF=alignments.lowResFluor2BF;
S2AHiRes=alignments.S2AHiRes;
Hi2LowResF=alignments.Hi2LowResF;
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;


%% set up low magvideos

aviFiles=dir([dataFolder filesep '20*.avi']);
aviFiles={aviFiles.name}';
aviFiles=aviFiles(cellfun(@(x) isempty(strfind(x,'HUDS')),aviFiles));
if length(aviFiles)==2
    aviFluorIdx=cellfun(@(x) ~isempty(strfind(x,'fluor')),aviFiles);
    %behaviorMovie=[dataFolder filesep aviFiles{~aviFluorIdx}];
    fluorMovie=[dataFolder filesep aviFiles{aviFluorIdx}];
else
    display('Select avi files, behavior and then low mag fluor');
    movies=uipickfiles('FilterSpec',dataFolder);
    %behaviorMovie=movies{1};
    fluorMovie=movies{2};
end

%behaviorVidObj = VideoReader(behaviorMovie);
fluorVidObj= VideoReader(fluorMovie);



%% set up high mag videos
rows=1200;
cols=600;
nPix=rows*cols;
if isempty(vidInfo)

[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);
else
    bfAll=vidInfo.bfAll;
    fluorAll=vidInfo.fluorAll;
    hiResData=vidInfo.hiResData;
    
end




%% set up timing alignments and lookups
bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');




%% load centerline
behaviorFolder=dir([dataFolder filesep 'Behavior*']);
behaviorFolder=behaviorFolder([behaviorFolder.isdir]);
behaviorFolder=[dataFolder filesep behaviorFolder(1).name];
centerlineFile=dir([behaviorFolder filesep 'center*']);
centerlineFile=[behaviorFolder filesep centerlineFile(1).name];
centerline=load(centerlineFile);
CLfieldNames=fieldnames(centerline);
CLfieldIdx=cellfun(@(x) ~isempty(strfind(x,'centerline')),CLfieldNames);
CLoffsetIdx=cellfun(@(x) ~isempty(strfind(x,'off')),CLfieldNames);
if any(CLoffsetIdx)
CLoffset=centerline.(CLfieldNames{CLoffsetIdx});


else
   CLoffset=0; 
end

%centerline=centerline.centerline;
centerline=centerline.(CLfieldNames{CLfieldIdx});


%% load images

         
         Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);

    hiResIdx=find(hiResData.stackIdx==iStack)+ zOffset;
    zRange=hiResData.Z(hiResIdx-zOffset);
    
            fluorIdx=round(fluorIdxLookup(hiResIdx));
         fluorIdxRange=[min(fluorIdx) max(fluorIdx)];
               fluorFrame=read(fluorVidObj,fluorIdxRange);
%             bfFrame = read(behaviorVidObj,bfIdxRange);
             fluorFrame=squeeze(fluorFrame(:,:,1,:));
%             bfFrame=squeeze(bfFrame(:,:,1,:));
% %                   bfFrame=imwarp(bfFrame,invert(lowResFluor2BF.t_concord),...
% %                 'OutputView',lowResFluor2BF.Rsegment);
% %               bfFrame=imwarp(bfFrame,Hi2LowResF.t_concord,...
% %                 'OutputView',Hi2LowResF.Rsegment);
% %                bfFrame=bfFrame((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
% 
             fluorFrame2=imwarp(fluorFrame,Hi2LowResF.t_concord,...
                 'OutputView',Hi2LowResF.Rsegment);
             if all(size(fluorFrame2(:,:,1))==[rows cols]);
          fluorFrame2=fluorFrame2((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
cropFlag=1;
             else 
                 cropFlag=0;
             end

    status=fseek(Fid,2*(hiResIdx(1))*nPix,-1);
    pixelValues=fread(Fid,nPix*(length(hiResIdx)),'uint16',0,'l');
    hiResImage=reshape(pixelValues,rows,cols,length(hiResIdx));
    
    %% crop and align hi mag images
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    segmentChannel=pedistalSubtract(segmentChannel);
%     activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
%     activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
%     activityChannel=pedistalSubtract(activityChannel);
%ctrlPoints=fiducialPoints{iStack};




    %% filter to find center Z (some old things here)
   worm2=segmentChannel;     
worm3=bpass3(worm2,.5,[20 20 3]);
worm3Smooth=smooth3(worm3,'box',[15,15,1]);


zSize=size(segmentChannel,3);

segmentChannel4=normalizeRange(worm3Smooth)>.05;
midZplot3=zeros(1,zSize);
for i=1:zSize
    [midZx,midZy]=find(bwmorph(segmentChannel4(:,:,i),'remove'));
    midZplotVal=max(pdist([midZx midZy]));
    if ~isempty(midZplotVal)
    midZplot3(i)=midZplotVal;
    end
end
midZplot3=normalizeRange(midZplot3(:));
[~,midZ1]=max(smooth(midZplot3,3));
%[~,midZ2]=max(smooth(midZplot,3));
midZ=(mean(midZ1));
if midZ>zSize/2
    midZ=floor(midZ);
else
    midZ=ceil(midZ);
end

    


%% try fix centerline alignemnt by looking at lowmag fluor and finding a correction offset
fluorFrame3=normalizeRange(double(fluorFrame2));
fluorFrame3=double((fluorFrame3>(graythresh(fluorFrame2(:)))));
%fluorFrame3(fluorFrame3<.1)=.1;

fluorProj=normalizeRange(sum(fluorFrame3,3));
hiResProj=normalizeRange(sum(segmentChannel,3));
%fluorProj(fluorProj>.5)=.5;
%fluorProj(fluorProj<.1)=.1;
fluorProj=normalizeRange(fluorProj);
fluorProj2=convnfft(fluorProj,Sfilter2,'same');
%fluorProj2=smooth2a(fluorProj,21,21);

corrIm=conv2(fluorProj,rot90(hiResProj,2),'same');
[CLoffsetY,CLoffsetX]=find(corrIm==max(corrIm(:)));
CLoffsetX=CLoffsetX-round(size(fluorProj,2)/2);
CLoffsetY=CLoffsetY-round(size(fluorProj,1)/2);

fluorProj2(fluorProj2<0)=0;
hiResIdxStretch=min(hiResIdx)-CLsearchWindow:max(hiResIdx)+CLsearchWindow;
bfIdx=round(bfIdxLookup(hiResIdxStretch));
CLIdx=bfIdx-CLoffset;
CL=centerline(:,:,CLIdx);


CLIdx=CLIdx-min(CLIdx)+1;
[~,ia,ic]=unique(CLIdx);

CL2=[];
[CL2(:,2,:),CL2(:,1,:)]=transformPointsInverse(lowResFluor2BF.t_concord,CL(:,2,:),CL(:,1,:));

%[CL2(:,1,:),CL2(:,2,:)]=transformPointsForward(Hi2LowResF.t_concord,CL2(:,2,:),CL2(:,1,:));
[CL2(:,1,:),CL2(:,2,:)]=transformPointsForward(Hi2LowResF.t_concord,CL2(:,2,:),CL2(:,1,:));
if cropFlag
CL2(:,2,:)=CL2(:,2,:)-(rect1(2)-1);
end
%[CL2(:,1,:),CL2(:,2,:)]=transformPointsForward(hiResFix.t_concord,CL2(:,2,:),CL2(:,1,:));

CL2X=reshape(CL2(:,1,CLsearchWindow:end-CLsearchWindow),[],1,1);
CL2Y=reshape(CL2(:,2,CLsearchWindow:end-CLsearchWindow),[],1,1);
[xyOffset2,~]=fminsearch(@(x) CLsearch(fluorProj2,CL2X+x(1),CL2Y+x(2),show),lastOffset);
%limit translation fix
% lastOffset=xyOffset3;
%   lastOffset(lastOffset>50)=50;
%   lastOffset(lastOffset<-50)=-50;

xyOffset3=xyOffset2-[CLoffsetX CLoffsetY];
  %xyOffset3(xyOffset3>50)=50;
  %xyOffset3(xyOffset3<-50)=-50;
CL2(:,2,:)=CL2(:,2,:)+xyOffset3(2);
CL2(:,1,:)=CL2(:,1,:)+xyOffset3(1);

if show
    close all
    imagesc(worm2(:,:,midZ))
    hold on
    CL2X=CL2X+xyOffset3(1);
    CL2Y=CL2Y+xyOffset3(2);
    
    plot(CL2X,CL2Y,'xr');
end




%% interpolate to parameterize by length

%
% CL3all=CL2(1:30,:,:);
% CL3all=interp1(CL3all,linspace(1,size(CL3all,1),300),'pchip');
CLlengthRange=2500;
CL3all=zeros(CLlengthRange,2,size(CL2,3));

%reinterpolate centerline by length


for iCL=1:size(CL2,3)
    CL3temp=CL2(:,:,iCL);
    s=[0; cumsum(squeeze(sqrt(sum((diff(CL3temp,[],1)).^2,2))))];
%CL3temp=interp1(s,CL3temp,0:10:40000);
CL3temp=interp1(s,CL3temp,1:1:CLlengthRange,'linear','extrap');
    CL3all(:,:,iCL)=CL3temp;
    if show
        if iCL==1
            close all
            imagesc(segmentChannel(:,:,midZ))
            hold on
        end
        plot(CL3temp(:,1),CL3temp(:,2))
    drawnow
    end
    
end


%% align centerlines parameterizations by correlation
% align centerlines using plotMatch (correlations) to account for
% centerline sliding
shiftVec=zeros(1,length(ia)-1);
for i=2:length(ia);
    CL1temp=CL3all(:,:,ia(i));
    CL2temp=CL3all(:,:,ia(i-1));
  %  CL1temp=bsxfun(@minus,CL1temp,mean(CL1temp));
   % CL2temp=bsxfun(@minus,CL2temp,mean(CL2temp));
    
    [corrtemp,r]=plotMatch(CL1temp,CL2temp,600);
    %[shift,~]=find(corrtemp==max(corrtemp(:)));
    [~,shift]=min(corrtemp);
    shiftVec(i)=r(shift);
%[i r(shift) shift];
end

shiftVec=cumsum(shiftVec);
shiftVec=shiftVec-shiftVec(round(length(shiftVec)/2));

CL3all2=zeros(CLlengthRange+501,2,size(CL2,3));

for iCL=1:size(CL2,3)
    CL3temp=CL3all(:,:,iCL);
CL3temp=interp1(CL3temp,shiftVec(ic(iCL))+(-500:CLlengthRange),'linear','extrap');
    CL3all2(:,:,iCL)=CL3temp;
    
    if show
        if iCL==1
            close all
            imagesc(segmentChannel(:,:,midZ))
            hold on
        end
        plot(CL3temp(:,1),CL3temp(:,2))
    drawnow
    end
    
    
end
 % center around head portion of worm that is in image
% inImage=sum(all(CL3all2>0 & CL3all2<600,2),3);
% inImage=inImage>max(inImage)/2;
% cc_inImage=bwconncomp(inImage);

CLcenter=sum(bsxfun(@minus, CL3all2,rect1(3:4)/2).^2,2);
[~,CLcenter]=min(CLcenter,[],1);


inImageRange=mean(CLcenter(:))+(-outputLength:outputLength);


%%
CL3all=interp1(CL3all2,inImageRange,'*linear','extrap');

%xyzs_avg = ActiveRevolutionFit(im, cline_para, CL3);

midIm=normalizeRange(mean(worm2(:,:,midZ+(-1:1)),3));
midIm=double(midIm);
midIm=bpass(midIm,2,[20,20]);
midImS=imfilter(midIm,Sfilter);
%midImS2=imfilter(midIm,Sfilter2);


%% more xy fixes
%minSearch=zeros(size(CL3all,1),length(ia));
CL3allX=(CL3all(:,1,CLsearchWindow:end-CLsearchWindow));
CL3allY=(CL3all(:,2,CLsearchWindow:end-CLsearchWindow));

minSearch=interp2(midImS,CL3allX,CL3allY);
minSearch=squeeze(nansum(minSearch,1));
minY=find(minSearch==max(minSearch));
midZCL=round(mean(minY))+CLsearchWindow;
CL3=CL3all(:,:,midZCL);


%%

if show
    close all
    
imagesc(midIm)

hold on
plot(CL3(:,1),CL3(:,2),'x')
axis equal
hold off
drawnow
end
% subplot(1,2,2)
% imagesc(midImST)
% hold on
% plot(CL3(:,1),CL3(:,2),'x')
% axis equal
% hold off


%% make coordinate system around the worm
pixel_interp=1;
Tv=zeros(size(CL3all,1),3,size(CL3all,3));
Bv=Tv; Nv=Tv;
for iSlice=1:size(CL3all,3);
% [Tv(:,:,iSlice),Bv(:,:,iSlice),Nv(:,:,iSlice)]=...
%     tbnVector([CL3all(:,1,iSlice),CL3all(:,2,iSlice),ones(size(CL3(:,1)))]);


T=normr(gradient(CL3all(:,:,iSlice)',5)');
N=[T(:,2) -T(:,1)];
B=T(:,1).*N(:,2)-T(:,2).*N(:,1);
N=bsxfun(@times, N,sign(B));
B=sign(B);

Tv(:,:,iSlice)=[T zeros(size(CL3(:,1)))];
Nv(:,:,iSlice)=[N zeros(size(CL3(:,1)))];
Bv(:,:,iSlice)=[zeros(size(CL3(:,1))) zeros(size(CL3(:,1))) B];


end
% 
% Bv=[Bv(1,:,:);Bv;];
% Nv=[Nv(1,:,:);Nv;];
signVector=sign(Bv(:,3,:));
Bv=bsxfun(@times,Bv,signVector);
Nv=bsxfun(@times,Nv,signVector);

%select worm orientation and fix 
if isempty(side)
    imagesc(midIm)
    choice = menu('Which Side is the nerve chord on?','Right','Left');
    if choice==2
        side='Left';
        Bv=-Bv;
        Nv=-Nv;
    else
        side='Right';
        
    end
else
    if ~isempty(strfind(side,'eft'))
          Bv=-Bv;
        Nv=-Nv;
    end
end




  plane_num=size(Tv,1);
        %make the first and last 'endround' tbn vectors the same so nothing
        %strange happens at the ends.
        if plane_num>10
            endround=5;
        else
            endround=round(plane_num/2);
        end
for iSlice=1:size(CL3all,3);
        for i=1:endround
            Tv(plane_num-i+1,:,iSlice)=Tv(plane_num-endround+1,:,iSlice);
            Bv(plane_num-i+1,:,iSlice)=Bv(plane_num-endround+1,:,iSlice);
            Nv(plane_num-i+1,:,iSlice)=Nv(plane_num-endround+1,:,iSlice);
            Tv(i,:,iSlice)=Tv(endround,:,iSlice);
            Bv(i,:,iSlice)=Bv(endround,:,iSlice);
            Nv(i,:,iSlice)=Nv(endround,:,iSlice);
        end
end

        %create a 2*window +1 square around each point for interpolationg using
        %the B and N vectors
        
        
  %% show stack with centerline
        if show
close all
% 
for iSlice=1:size(worm2,3);
    
    imagesc(worm2(:,:,iSlice));colormap hot
hold on
    clSlice=CLsearchWindow+iSlice;
    clSlice=round(clSlice);
    clSlice(clSlice<1)=1;
    clSlice(clSlice>size(CL3all,3))=size(CL3all,3);
%scatter(CL3all(:,1,clSlice),CL3all(:,2,clSlice),[],1:length(CL3all(:,1,clSlice)),'.');
plot(CL3all(:,1,clSlice),CL3all(:,2,clSlice));
 quiver(CL3all(1:10:end,1,clSlice),CL3all(1:10:end,2,clSlice),...
     Nv(1:10:end,1,clSlice),Nv(1:10:end,2,clSlice))

hold off
axis auto equal
xlim([0 600]);ylim([0 600])
%print(gcf,['Y:\Jeff\PowerPoint\New folder\MySavedPlot' num2str(iSlice,'%3.5d') ],'-dpng')
pause(.1)
end
        end

 %% straighten interpolation

        [J,K]=meshgrid(-outputRadius:1/pixel_interp:outputRadius,...
            -outputRadiusZ:1/pixel_interp:outputRadiusZ);


        zslice=bsxfun(@times,J,permute(Nv(:,3,1),[3,2,1]))*zRatio+...
            bsxfun(@times,K,permute(Bv(:,3,1),[3,2,1]))*zRatio+midZ;
        
        
zLevels=((zslice(:,1,1)));

%correct for non monoticity
if sign(nanmean(diff(zRange)))==-1
    
    zRange2=unique(cummin(zRange));
    zInterp=interp1(zRange2-zRange2(midZ),1:length(zRange2),zLevels/10-zLevels(round(outputRadiusZ+1))/10);
zInterp=flipud(zInterp);
  zLevels=flipud(zLevels);  
else
    zRange2=unique(cummax(zRange));

    zInterp=interp1(zRange2-zRange2(midZ),1:length(zRange2),zLevels/10-zLevels(round(outputRadiusZ+1))/10);

end


zslice=repmat(zInterp,1,2*outputRadius+1,size(Bv,1));

zLevels(zLevels<min(ia))=min(ia);
zLevels(zLevels>size(worm2,3))=size(worm2,3);
zLevels=zLevels+CLsearchWindow;
CLzLevels=zLevels;
% zLevels=interp1(ia,CLIdxUnique,zLevels,'linear','extrap');
% zLevels(zLevels<1)=1;
% zLevels(zLevels>max(ia))=max(ia);
CL3xinterp=interp1(squeeze(CL3all(:,1,:))',CLzLevels,'linear')';
CL3xinterp=permute(CL3xinterp,[2,3,1]);
CL3yinterp=interp1(squeeze(CL3all(:,2,:))',CLzLevels,'linear')';
CL3yinterp=permute(CL3yinterp,[2,3,1]);
NvInterpx=interp1(squeeze(Nv(:,1,:))',CLzLevels,'linear');
NvInterpy=interp1(squeeze(Nv(:,2,:))',CLzLevels,'linear');
BvInterpx=interp1(squeeze(Bv(:,1,:))',CLzLevels,'linear');
BvInterpy=interp1(squeeze(Bv(:,2,:))',CLzLevels,'linear');

        xslice=bsxfun(@times,J,permute(NvInterpx,[1,3,2]))+...
            bsxfun(@times,K,permute(BvInterpx,[1,3,2]));
        xslice=bsxfun(@plus, xslice,CL3xinterp);
         yslice=bsxfun(@times,J,permute(NvInterpy,[1,3,2]))+...
            bsxfun(@times,K,permute(BvInterpy,[1,3,2]));
        yslice=bsxfun(@plus, yslice,CL3yinterp);
        
xslice=permute(xslice,[3,2,1]);
yslice=permute(yslice,[3,2,1]);
zslice=permute(zslice,[3,2,1]);

%%
        %use points to interpolate, XY in matlab is messed me up.. but this
        %works
        %if using gradient, convolve stack with sobel operator in each
        %direction and find magnitude
        
 %       if 1%~cline_para.gradflag
 
 xslice=round(xslice);yslice=round(yslice);zslice=round(zslice);
 inImageMap=xslice>0 & zslice>0 & yslice>0 & xslice<size(worm2,2) &...
     yslice<size(worm2,1) &  zslice<size(worm2,3);
 inImageMapIdx=sub2ind_nocheck(size(worm2),(yslice(inImageMap)),...
     (xslice(inImageMap)),(zslice(inImageMap)));
 V=zeros(size(xslice));
 V(inImageMap)=worm2(inImageMapIdx);

     %       Vsmooth=interp3(worm3,xslice,yslice,zslice,'*linear',0);
            
 %           A=interp3(activityChannel,xslice,yslice,zslice,'*nearest',0);
%         else          lowResV=interp3(fluorFrame2,xslice,yslice,zslice,'*linear',0);

   %% stack stabilization
   
   [~, tformAll]=stackStabilization(V,30,show,0);
R = imref2d(size(V(:,:,1))) ;
   for iSlice=1:size(xslice,3)
       if any(any(inImageMap(:,:,iSlice)))
temp=imwarp(cat(3,xslice(:,:,iSlice),yslice(:,:,iSlice)),R,tformAll{iSlice},'nearest',...
                'OutputView',R);
            xslice(:,:,iSlice)=temp(:,:,1);
            yslice(:,:,iSlice)=temp(:,:,2);

       end
   end

 xslice=round(xslice);yslice=round(yslice);zslice=round(zslice);
 inImageMap=xslice>0 & zslice>0 & yslice>0 & xslice<size(worm3,2) &...
     yslice<size(worm3,1) &  zslice<size(worm3,3);
 inImageMapIdx=sub2ind(size(worm3),(yslice(inImageMap)),...
     (xslice(inImageMap)),(zslice(inImageMap)));
 Vsmooth=zeros(size(xslice));
  V=Vsmooth;
 Vsmooth(inImageMap)=worm3(inImageMapIdx);
  V(inImageMap)=worm2(inImageMapIdx);

   %%
   

Vproj=squeeze(nansum(V,3));


%% Correlation algin with template image

if ~isempty(Vtemplate)
% repalced xcorr with conv2 for small speed boost, search area is decreased
xIm=conv2(Vproj,rot90(Vtemplate,2),'same');
[xlag,ylag]=find(xIm==max(xIm(:)));
lags=[xlag,ylag]-round(size(Vproj)/2);

[ndX,ndY,ndZ]=ndgrid(1:size(V,1),1:size(V,2),1:size(V,3));
ndX=ndX+lags(1);
ndY=ndY+lags(2);
ndX=round(ndX);ndY=round(ndY);
inImage=(ndY>0 & ndX>0 & ndX<(2*outputLength+1) & ndY<(2*outputRadius+1));
inImageIdx=sub2ind_nocheck(size(V),ndX(inImage),ndY(inImage),ndZ(inImage));
temp=zeros(size(V));
temp(inImage)=V(inImageIdx);
V=temp;
temp(inImage)=xslice(inImageIdx);
xslice=temp;
temp(inImage)=yslice(inImageIdx);
yslice=temp;
temp(inImage)=zslice(inImageIdx);
zslice=temp;
temp(inImage)=Vsmooth(inImageIdx);
Vsmooth=temp;

else
Vproj=squeeze(nansum(V,3));
end




%% now segmentation

    V(isnan(V))=0;
    Vsmooth(isnan(Vsmooth))=0; % option, to use presmoothed version, much faster but may not be a s good
    imsize=size(V);
    [wormBW2,~]=WormSegmentHessian3dStraighten(V,options,Vsmooth);
     %%
    BWplot=(squeeze(sum(sum(wormBW2,1),2)));
    BWplot=smooth(BWplot,20);
    [~,locs]=findpeaks(BWplot);
    endpts=locs([1,length(locs)]);
    [~,locs]=findpeaks(-BWplot);
    botpoint1=locs((locs>endpts(1)));
    if isempty(botpoint1);botpoint1=1;end;
    botpoint1=botpoint1(1);
    botpoint2=locs((locs<endpts(2)));
    if isempty(botpoint2);botpoint2=imsize(3);end;
    
    botpoint2=botpoint2(end);
    botpoint1(botpoint1>imsize(3)*1/4)=1;
    botpoint2(botpoint2<imsize(3)*3/4)=imsize(3);
    
    
    cc=bwconncomp(wormBW2,6);
    
    badRegions=(cellfun(@(x) any(zindexer(x,imsize(1)*imsize(2))<=botpoint1),cc.PixelIdxList)...
        |cellfun(@(x) any(zindexer(x,imsize(1)*imsize(2))>=botpoint2),cc.PixelIdxList))';
    
    wormBW2(cell2mat(cc.PixelIdxList(badRegions)'))=false;
    cc.PixelIdxList=cc.PixelIdxList(~badRegions);
    cc.NumObjects=nnz(~badRegions);
     cc=bwconncomp(wormBW2,6);
    stats=regionprops(cc,V,'Centroid','MeanIntensity',...
        'Area');
    
     intensities=[stats.MeanIntensity]';  
    P=[cell2mat({stats.Centroid}'),iStack*ones(cc.NumObjects,1)...
        (1:cc.NumObjects)'  intensities];
    P(:,[1 2])=P(:,[2 1]);
Areas=[stats.Area]';
    Poriginal=[interp3(xslice,P(:,2),P(:,1),P(:,3)) ...
        interp3(yslice,P(:,2),P(:,1),P(:,3)) ...
        interp3(zslice,P(:,2),P(:,1),P(:,3))];
    pointStats.straightPoints=P(:,1:3);
    pointStats.rawPoints=Poriginal;
    pointStats.stackIdx=iStack;
    pointStats.pointIdx=(1:cc.NumObjects)';
    pointStats.Rintensities=intensities;
    pointStats.Volume=Areas;
    pointStats.baseImg=logical(wormBW2);
    pointStats.transformx=uint16(xslice);
    pointStats.transformy=uint16(yslice);
    pointStats.transformz=uint16(zslice);
 %   pointStats(counter).Gintensities=intensities;
 
 %% tran
if ~isempty(ctrlPoints)
ctrlPoints2=nan(size(ctrlPoints,1),3);
fiducialIdx=cell2mat(cellfun(@(x) ~isempty(x),ctrlPoints(:,1),'uniform',0));

ctrlPoints2(fiducialIdx,:)=cell2mat(ctrlPoints(:,[1 2 4]));
ctrlPoints2(:,3)=ctrlPoints2(:,3)-min(hiResIdx)+1;
pointStats.controlPoints=ctrlPoints2;
    
else
pointStats.controlPoints=[];
end

    
    
    %%    
    
    %TrackData{counter}=P;
    

%fileName=[imageFolder filesep 'image' num2str(iStack,'%3.5d') '.tif'];
fileName2=[imageFolder2 filesep 'image' num2str(iStack,'%3.5d') '.tif'];
fileName3=[imageFolder2 filesep 'pointStats' num2str(iStack,'%3.5d')];

%tiffwrite(fileName,Vproj,'tif');
tiffwrite(fileName2,single(V),'tif');
save(fileName3,'pointStats'); 
%save with compression reduces file size by more than 70%

%save(fileName3,'wormRegions');

fclose(Fid);
catch me
    fileName=[imageFolder2 filesep 'ERROR' num2str(iStack,'%3.5d')];
    save(fileName);
    rethrow(me)
end