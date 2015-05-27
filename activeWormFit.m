% worm=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResSegmentFolder3Dtest_raw\image00500.tif');
% activity=stackLoad('E:\20141212\BrainScanner20141212_145951\fiducials\hiResActivityFolder3Dtest_raw\image00500.tif');

% dataFolder=uipickfiles;
% dataFolder=dataFolder{1};
dataFolder='F:\20141212\BrainScanner20141212_145951\';
%[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,imSize);
%%





rows=1200;
cols=600;
nPix=rows*cols;
if exist([dataFolder filesep 'hiResData.mat'],'file')
    hiResData=load([dataFolder filesep 'hiResData']);
    hiResData=hiResData.dataAll;
else
    hiResData=highResTimeTraceAnalysisTriangle4(dataFolder,rows,cols);
end

[bfAll,fluorAll,hiResData]=tripleFlashAlign(dataFolder,[rows cols]);

bfIdxList=1:length(bfAll.frameTime);
fluorIdxList=1:length(fluorAll.frameTime);
bfIdxLookup=interp1(bfAll.frameTime,bfIdxList,hiResData.frameTime,'linear');
fluorIdxLookup=interp1(fluorAll.frameTime,fluorIdxList,hiResData.frameTime,'linear');

[hiImageIdx,ib]=unique(hiResData.imageIdx);
hiResLookup=interp1(hiImageIdx,ib,1:length(hiResData.frameTime));



%% load alignment data

display('Select Low Res Alignment')
%lowResFluor2BF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
lowResFluor2BF=load('Y:\CommunalCode\3dbrain\registration\20141212LowResBehavior2Fluor.mat');
lowResBF2FluorT=invert(lowResFluor2BF.t_concord);


display('Select Hi to Low Fluor Res Alignment')
%Hi2LowResF=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
Hi2LowResF=load('Y:\CommunalCode\3dbrain\registration\20141212HighResS2LowResFluorBeads.mat');


% display('Select Hi to Low Res Alignment')
% 
% Hi2LowRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
% Hi2LowRes=load(Hi2LowRes{1});
% t_concord = fitgeotrans(Hi2LowRes.Sall,Hi2LowRes.Aall,'projective');
 display('Select Hi Res Alignment')

%S2AHiRes=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\registration');
S2AHiRes=load('Y:\CommunalCode\3dbrain\registration\20141212HiResS2A.mat');
rect1=S2AHiRes.rect1;
rect2=S2AHiRes.rect2;


%% load Fiducials file
fiducialFile=dir([dataFolder filesep '*iducial*']);
fiducialFile={fiducialFile.name}';
if length(fiducialFile)~=1
        display('Select model file');

    fiducialFile=uipickfiles('FilterSpec',dataFolder);
    fiducialFile=load(fiducialFile{1});
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;
    
else
    fiducialFile=load([dataFolder filesep fiducialFile{1}]);
    fiducialPoints=fiducialFile.fiducialPoints;
    z2ImageIdxOffset=fiducialFile.timeOffset;

end

%%

gaussKernal=gausswin(30);
gaussKernal=convnfft(gaussKernal,gaussKernal');
gaussKernal=convnfft(gaussKernal,permute(gausswin(3),[2,3,1]));
gaussKernal=gaussKernal/sum(gaussKernal(:));
%%
gaussKernal2=gausswin(200);
gaussKernal2=convnfft(gaussKernal2,gaussKernal2');
%%
gaussKernal50=gausswin(30);
gaussKernal50=convnfft(gaussKernal50,gaussKernal50');
%gaussKernal50=del2(gaussKernal50);
gaussKernal50=gradient(gaussKernal50)';
gaussKernal50(gaussKernal50<0)=gaussKernal50(gaussKernal50<0)*10;
%gaussKernal50(gaussKernal50>0)=gaussKernal50(gaussKernal50>0)*3;


% gaussx=meshgrid(1:length(gaussKernal50))';
% gaussx=gaussx-mean(gaussx(:));
% gaussKernal50(gaussKernal50>0 & gaussx>0)=0;

%%

sphereKernal=gausswin(80);
sphereKernal=convnfft(sphereKernal,sphereKernal');
sphereKernal=convnfft(sphereKernal,permute(gausswin(20),[2,3,1]));
sphereKernal=sphereKernal/max(sphereKernal(:));
sphereKernal(sphereKernal>.5)=0;
sphereKernal(sphereKernal>.2)=1;

%%
cline_para.CLalpha=10;
cline_para.CLbeta=10;
cline_para.kappa=5;
cline_para.CLalphaR=100;
cline_para.CLbetaR=1;
cline_para.kappaR=1;
cline_para.gamma=10;
cline_para.zRatio=1/3;
cline_para.gradient_force=1;
cline_para.iterations=300;
cline_para.show_flag=2;


%%%
% CL3=[CL2(:,1)+50,CL2(:,2),0*CL2(:,1)+1,0*CL2(:,1)];
% CL3=CL3(1:20,:);
% imZproj=sum(im,3);
% 
% xyzs_avg = ActiveRevolutionFit(imZproj, cline_para, CL3);


        %%
        Options.Registration='Affine';
  Options.Similarity='p';
  Options.MaxRef=20;
Options.Verbose=0;
Options.SigmaFluid=6;

%%
centerline=load([dataFolder filesep 'Behavior Analysis' filesep 'centerline']);
CLoffset=centerline.offset;
centerline=centerline.centerline;


%%
zOffset=z2ImageIdxOffset;

startStack=500;
endStack=800;
se=strel('disk',30);
imageFolder=[dataFolder filesep 'straightProj2'];
imageFolder2=[dataFolder filesep 'CLstraight4'];
mkdir(imageFolder);
mkdir(imageFolder2);
counter=1;

%%
% parfor_progress(0)
% parfor_progress(endStack-startStack);
stackRange= 499+  [163  164];

progressbar
for iStack=stackRange
    counter=1+iStack-startStack;
progressbar((iStack-startStack)/(endStack-startStack));
 tic
     try
         Fid=fopen([dataFolder filesep 'sCMOS_Frames_U16_1024x1024.dat']);


   segmentChannel2=[]; segmentChannel3=[]; segmentChannel4=[]; segmentChannel5=[];

    hiResIdx=find(hiResData.stackIdx==iStack)-zOffset;
  %  progressbar(iStack/length(annotated));

  %    currentFiducials=fiducialPoints{annotated(iStack)};
% 
%         fiducialIdx=find(cell2mat((cellfun(@(x)( ~isempty(x)),currentFiducials(:,1),'uniformoutput',0))));
%         fiducialIdx(isnan(cell2mat(currentFiducials(fiducialIdx,1))))=[];
%           fiducialZ=round(cell2mat(currentFiducials(fiducialIdx,4)));
% 
%         currentFiducials=cell2mat(currentFiducials(fiducialIdx,:));    
%         Rout=zeros(12,1);
%         Gout=Rout;
%         backgroundOut=Gout;
%         
%         rFiducialPoints=round(currentFiducials(:,[1,2,4]));
  
%%
    
%progressbar(iStack/length(annotated),i/size(rFiducialPoints,1));


%    inPlaneFiducials=currentFiducials(ib==i,:);
    
    status=fseek(Fid,2*(hiResIdx(1))*nPix,-1);
    pixelValues=fread(Fid,nPix*(length(hiResIdx)),'uint16',0,'l');
 
    
    hiResImage=reshape(pixelValues,rows,cols,length(hiResIdx));
    
    segmentChannel=hiResImage((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3),:);
    activityChannel=hiResImage((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3),:);
    activityChannel=imwarp(activityChannel,S2AHiRes.t_concord,'OutputView',S2AHiRes.Rsegment);
    segmentChannel=pedistalSubtract(segmentChannel);
    activityChannel=pedistalSubtract(activityChannel);
    
ctrlPoints=fiducialPoints{iStack};
    %%
   worm2=segmentChannel;     
worm3=bpass3(worm2,.3,[50 50 3]);

segmentChannel3=convnfft(worm3,sphereKernal,'same');

sumIm=normalizeRange(nanmean(worm2,3));
sumIm=smooth2a(sumIm,15,15);
sumIm=normalizeRange(sumIm);

%%
zSize=size(segmentChannel,3);

for i=1:zSize;
    segmentChannel2(:,:,i)=imtophat(-smooth2a(worm3(:,:,i),15,15),se);

end
segmentChannel2=normalizeRange(segmentChannel2);


%segmentChannel4=bsxfun(@times,segmentChannel2,sumIm);


%%
segmentChannel4=(segmentChannel2>graythresh(segmentChannel2(:)));
%%

for i=1:zSize;
    segmentChannel5(:,:,i)=bwdist(~segmentChannel4(:,:,i));

end

%midZplot=double(squeeze(max(max(segmentChannel5,[],1),[],2)));
midZplot=double(squeeze(sum(sum(segmentChannel5>max(segmentChannel5(:))/2,1),2)));
%%
%midZplot=squeeze(mean(nanmean(segmentChannel4,1),2));
midZplot3=squeeze(mean(nanmean(segmentChannel3,1),2));
midZplot3=normalizeRange(midZplot3);
pl=smooth(midZplot,1);
pl=normalizeRange(pl);
[pks,locs,w,p] = findpeaks(pl);
if any(pks>.75)
locs(pks<.75)=[];
pks(pks<.75)=[];
end
%pks=midZplot3(locs);
if any(~(locs<10 |locs>30))
pks(locs<10 |locs>30)=0;
end


%midZ=locs(pks==max(pks));

    [~,locPos]=max(pks);
    midZ=locs(locPos);

if isempty(midZ)
    midZ=round(size(segmentChannel,3)/2);
end



%%
bfIdx=round(bfIdxLookup(hiResIdx));
CLIdx=bfIdx-CLoffset;
CL=centerline(:,:,CLIdx);


CLIdx=CLIdx-min(CLIdx)+1;
ia2=accumarray(CLIdx,1:length(CLIdx),[],@mean);
[CLIdxUnique,ia,ic]=unique(CLIdx);

CL2=[];
[CL2(:,2,:),CL2(:,1,:)]=transformPointsInverse(lowResFluor2BF.t_concord,CL(:,2,:),CL(:,1,:));

[CL2(:,1,:),CL2(:,2,:)]=transformPointsForward(Hi2LowResF.t_concord,CL2(:,2,:),CL2(:,1,:));
 CL2(:,2,:)=CL2(:,2,:)-600;


%%

%
% CL3all=CL2(1:30,:,:);
% CL3all=interp1(CL3all,linspace(1,size(CL3all,1),300),'pchip');
CL3all=[];

%reinterpolate centerline by length

for iCL=1:size(CL2,3)
    CL3temp=CL2(:,:,iCL);
    s=[0; cumsum(squeeze(sum((diff(CL3temp,[],1)).^2,2)))];
CL3temp=interp1(s,CL3temp,133.333:133.333:40000);
    CL3all(:,:,iCL)=CL3temp;
    
    
    
end

%align centerlines using plotMatch
shiftVec=0;
for i=2:length(ia);
    CL1temp=CL3all(:,:,ia(i));
    CL2temp=CL3all(:,:,ia(i-1));
  %  CL1temp=bsxfun(@minus,CL1temp,mean(CL1temp));
   % CL2temp=bsxfun(@minus,CL2temp,mean(CL2temp));
    
    [corrtemp,r]=plotMatch(CL1temp,CL2temp);
    %[shift,~]=find(corrtemp==max(corrtemp(:)));
    [~,shift]=min(corrtemp);
    shiftVec(i)=r(shift);
end
shiftVec=cumsum(shiftVec);
%%
CL3all2=[];

for iCL=1:size(CL2,3)
    CL3temp=CL3all(:,:,iCL);
CL3temp=interp1(CL3temp,shiftVec(ic(iCL))+(1:300),'linear','extrap');
    CL3all2(:,:,iCL)=CL3temp;
    
    
    
end

CL3all=CL3all2;
CL3=CL3all(:,:,midZ);

    
    

%%
Sfilter=max(gaussKernal2(:))-gaussKernal2;
 %Sfilter(Sfilter<.1)=-1;
% Sfilter(Sfilter>.8)=nan;
% Sfilter(isnan(Sfilter))=0;

%Sfilter(Sfilter>0)=nnz(Sfilter<0)/nnz(Sfilter>0);
% 
Sfilter(Sfilter<.1)=-(.1-Sfilter(Sfilter<.1))*80;%Sfilter(Sfilter<.01)-.3;

Sfilter(Sfilter>.8)=0;
Sfilter(Sfilter>0)=1;%nnz(Sfilter<0)/nnz(Sfilter>0);
%Sfilter=Sfilter-mean(Sfilter(:));

Sfilter2=max(gaussKernal2(:))-gaussKernal2;
Sfilter2(Sfilter2<.01)=Sfilter2(Sfilter2<.01)-.3;
Sfilter2(Sfilter2>.6)=0;

Sfilter2(Sfilter2>0)=nnz(Sfilter2<0)/nnz(Sfilter2>0);
Sfilter2=Sfilter2-mean(Sfilter2(:));

%xyzs_avg = ActiveRevolutionFit(im, cline_para, CL3);

midIm=normalizeRange(mean(worm2(:,:,midZ+[-1:1]),3));
midIm=double(midIm);
midIm=bpass(midIm,2,[20,20]);
midImS=imfilter(midIm,Sfilter);
midImS2=imfilter(midIm,Sfilter2);
%midImS(imdilate(sumIm<.09,true(10)))=0;%min(midImS(:))/2;


midImST=imtophat(-midIm,se);
midImST=normalizeRange(midImST);

midImST=midImST>.1;
midImST=imerode(midImST,true(10));
cc=bwconncomp(midImST);
areas=cell2mat(cellfun(@(x) length(x),cc.PixelIdxList,'uniform',0));
maxArea=find(areas==max(areas),1,'first');
midImSTBW=false(size(midImST));
midImSTBW(cc.PixelIdxList{maxArea})=true;
midImSTBW=bwmorph(midImSTBW,'shrink',Inf);
stats=regionprops(midImSTBW,'Centroid');
offset2=round(stats.Centroid);

%[offset2(2) offset2(1)]=find((midImST.*midImS)==max(midImST(:).*midImS(:)));
%xyOffset=fminsearch(@(x) CLsearch(midImST,CL3(:,1)+x(1),CL3(:,2)+x(2)),[0 0]);
% CL3(:,1)=CL3(:,1)+xyOffset(1);
% CL3(:,2)=CL3(:,2)+xyOffset(2);
%%
minSearch=[];
for i=1:length(CL3)
    minSearch(i)=CLsearch(midImS,CL3(:,1)-CL3(i,1)+offset2(1),...
        CL3(:,2)-CL3(i,2)+offset2(2),show);
    
%         minSearch(i)=CLsearch(sum(segmentChannel2,3),CL3(:,1)-CL3(i,1)+offset2(1),...
%         CL3(:,2)-CL3(i,2)+offset2(2),1);

end
%%
%minSearch(1:100)=Inf;
minSearch=smooth(minSearch,4);
[pks,locs]=findpeaks(-minSearch);
pks=-pks;
%[~,minIdx]=min(minSearch);
xyOffset=bsxfun(@plus,-CL3(locs,1:2),offset2);
v=sqrt(sum(xyOffset.^2,2));
if any (v<400)
minIdx=find(v<400);
pks=pks(minIdx);

minIdx=minIdx(pks==min(pks));
else 
   [~, minIdx]=min(v);
end
%minIdx=locs(minIdx);
%
xyOffset=xyOffset(minIdx,:);
[xyOffset3,fval]=fminsearch(@(x) CLsearch(midImS,CL3(:,1)+x(1),CL3(:,2)+x(2),show),xyOffset);

% if fval>-10
%   xyOffset3=fmincon(@(x) CLsearch(midImS,CL3(:,1)+x(1),CL3(:,2)+x(2),1),xyOffset,...
%       [],[],[],[],xyOffset-60,xyOffset+60);
% end

%%

CL3(:,1)=CL3(:,1)+xyOffset3(1);
CL3(:,2)=CL3(:,2)+xyOffset3(2);
CL3all(:,1,:)=CL3all(:,1,:)+xyOffset3(1);
CL3all(:,2,:)=CL3all(:,2,:)+xyOffset3(2);
%offsetAll(iStack-startStack+1,:)=xyOffset3;
%%
minSearch2=[];
for i=1:size(CL3all,3)
    minSearch2(i)=CLsearch(segmentChannel3(:,:,1),CL3all(:,1,i),...
        CL3all(:,2,i),show);
    
%         minSearch(i)=CLsearch(sum(segmentChannel2,3),CL3(:,1)-CL3(i,1)+offset2(1),...
%         CL3(:,2)-CL3(i,2)+offset2(2),1);
end
minSearch2=smooth(minSearch2,4);

minZ=find(minSearch2==min(minSearch2(1:min(midZ-4,10))),1,'first');

for i=1:size(CL3all,3)
    minSearch2(i)=CLsearch(segmentChannel3(:,:,end),CL3all(:,1,i),...
        CL3all(:,2,i),show);
    
%         minSearch(i)=CLsearch(sum(segmentChannel2,3),CL3(:,1)-CL3(i,1)+offset2(1),...
%         CL3(:,2)-CL3(i,2)+offset2(2),1);
end
minSearch2=smooth(minSearch2,4);
maxZ=find(minSearch2==min(minSearch2(max(midZ+4,size(worm3,3)-10):end)),1,'first');




% clear temp;
% for i=1:zSize
% temp(:,i)=interp2(worm2(:,:,i),CL3(:,1),CL3(:,2));
% end
% segmentChannevlSTD=nansum(temp,1);
% segmentChannelSTD=detrend(segmentChannelSTD);
% segmentChannelSTD=segmentChannelSTD-min(segmentChannelSTD);
% segmentChannelSTD=normalizeRange(segmentChannelSTD);
% f=fit((1:zSize)',segmentChannelSTD','gauss2','StartPoint',[1,15,5,1,25,5]);
% 
%  botPlane=[min(f.b1,f.b2),max(f.b1,f.b2)];
%  
%  botPlane=round(botPlane);
%  [~,midZ]=min(f((botPlane(1)):botPlane(2)));
%  midZ=midZ+botPlane(1)-1

% midZ=round(botPlane)+[10,-10];
% if segmentChannelSTD(midZ(1))<segmentChannelSTD(midZ(2));
%     midZ=midZ(2);
% else
%     midZ=midZ(1);
% end

%midZ=(f.b1/f.c1+f.b2/f.c2)/(1/f.c1+1/f.c2);
%%
% midZ=round(midZ);
% midIm=normalizeRange(worm2(:,:,midZ));
% CL3(:,3)=midZ;
% %
% subplot(1,2,1)
if show
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


%%
window=100;
pixel_interp=1;
Tv=[]; Bv=[]; Nv=[];
for iSlice=1:size(CL3all,3);
[Tv(:,:,iSlice),Bv(:,:,iSlice),Nv(:,:,iSlice)]=...
    tbnVector([CL3all(:,1,iSlice),CL3all(:,2,iSlice),ones(size(CL3(:,1)))]);
end

Bv=[Bv(1,:,:);Bv;];
Nv=[Nv(1,:,:);Nv;];
signVector=sign(Bv(:,3,:));
Bv=bsxfun(@times,Bv,signVector);
Nv=bsxfun(@times,Nv,signVector);
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
        
        
        %%
        if show
close all
% 
for iSlice=1:size(worm2,3);
    
    imagesc(worm2(:,:,iSlice));
    hold on
    iSlice=interp1([1 midZ size(worm2,3)],[minZ midZ maxZ], iSlice);
    iSlice=round(iSlice);
    iSlice(iSlice<1)=1;
    iSlice(iSlice>size(worm2,3))=size(worm2,3);
scatter(CL3all(:,1,iSlice),CL3all(:,2,iSlice),[],1:length(CL3all(:,1,iSlice)),'.');
quiver(CL3all(1:10:end,1,iSlice),CL3all(1:10:end,2,iSlice),...
    Nv(1:10:end,1,iSlice),Nv(1:10:end,2,iSlice))

hold off
axis auto equal
pause(.1)
end
        end

        %%
        xslice=zeros(pixel_interp*2*window+1,...
            pixel_interp*2*window+1,size(Bv,1));
        yslice=xslice;
        zslice=xslice;
        
        [J,K]=meshgrid(-window:1/pixel_interp:window,...
            -window:1/pixel_interp:window);
        %%

        zslice=bsxfun(@times,J,permute(Nv(:,3,1),[3,2,1]))*cline_para.zRatio+...
            bsxfun(@times,K,permute(Bv(:,3,1),[3,2,1]))*cline_para.zRatio+midZ;
        
        
zLevels=(squeeze(zslice(:,1,1)));
zLevels(zLevels<min(ia))=min(ia);
zLevels(zLevels>max(ia))=max(ia);
zLevels=interp1([1 midZ max(zLevels)],[minZ midZ maxZ], zLevels);
zLevels=interp1(ia,CLIdxUnique,zLevels,'linear','extrap');
zLevels(zLevels<1)=1;
zLevels(zLevels>size(ia,1))=size(ia,1);
CL3xinterp=interp1(squeeze(CL3all(:,1,ia))',zLevels,'linear')';
CL3xinterp=permute(CL3xinterp,[2,3,1]);
CL3yinterp=interp1(squeeze(CL3all(:,2,ia))',zLevels,'linear')';
CL3yinterp=permute(CL3yinterp,[2,3,1]);
NvInterpx=interp1(squeeze(Nv(:,1,ia))',zLevels,'linear');
NvInterpy=interp1(squeeze(Nv(:,2,ia))',zLevels,'linear');
BvInterpx=interp1(squeeze(Bv(:,1,ia))',zLevels,'linear');
BvInterpy=interp1(squeeze(Bv(:,2,ia))',zLevels,'linear');

        xslice=bsxfun(@times,J,permute(NvInterpx,[3,1,2]))+...
            bsxfun(@times,K,permute(BvInterpx,[3,1,2]));
        xslice=bsxfun(@plus, xslice,CL3xinterp);
         yslice=bsxfun(@times,J,permute(NvInterpy,[3,1,2]))+...
            bsxfun(@times,K,permute(BvInterpy,[3,1,2]));
        yslice=bsxfun(@plus, yslice,CL3yinterp);
        
xslice=permute(xslice,[3,2,1]);
yslice=permute(yslice,[3,2,1]);
zslice=permute(zslice,[3,2,1]);

        
        %use points to interpolate, XY in matlab is messed me up.. but this
        %works
        %if using gradient, convolve stack with sobel operator in each
        %direction and find magnitude
        
 %       if 1%~cline_para.gradflag
           
            V=interp3(worm2,xslice,yslice,zslice,'*linear',0);
            A=interp3(activityChannel,xslice,yslice,zslice,'*linear',0);
%         else
%             kernal=fspecial('sobel');
%             kernal3(:,:,1)=kernal;
%             kernal3(:,:,2)=2*kernal;
%             kernal3(:,:,3)=kernal;
%             kernal3r=kernal3;
%             kernal3p=permute(kernal3,[2,3,1]);
%             kernal3c=permute(kernal3,[3,1,2]);
%             hostack=convn(im,kernal3r,'same');
%             vertstack=convn(im,kernal3c,'same');
%             pagestack=convn(im,kernal3p,'same');
%             pagestack(:,:,1)=0;
%             pagestack(:,:,end)=0;
%             
%             grad_stack=sqrt(hostack.^2+vertstack.^2+pagestack.^2);
%             V=interp3(grad_stack,yslice,xslice,zslice,'linear',0);
%             V(1,:,:)=0;
%             V(end,:,:)=0;
%             V(:,1,:)=0;
%             V(:,end,:)=0;
%             
%             V=smooth3(V,'box',[5,5,5]);
%        end
%         if flag.zplanar==1 && flag.straight==1
%             % V=smooth3(V,'box',[1,1,3]);
%             %        V=smooth3(V,'gaussian',[5,5,5],3*pixel_interp);
%         end
%         
%         if flag.showfits~=1
%             clear xlice yslice zslice
%         end




%% Unwrap, not currently used
% [t_vector,b_vector,n_vector]=tbnVector([CL3(:,1),CL3(:,2),0*CL3(:,1)]);
% b_vector=[b_vector(1,:);b_vector;];
% n_vector=[n_vector(1,:);n_vector;];
% t_vector=[t_vector(1,:);t_vector;];
% 
% signVector=sign(n_vector(:,1));
% %b_vector=bsxfun(@times,b_vector,signVector);
% n_vector=bsxfun(@times,n_vector,signVector);
% %%
% radii=1:140;
% radii=reshape(radii,1,1,[]);
% dtheta=pi/40;
% theta=dtheta:dtheta:2*pi;
% nx=(b_vector(:,1))*cos(theta)+(n_vector(:,1))*sin(theta);
% ny=(b_vector(:,2))*cos(theta)+(n_vector(:,2))*sin(theta);
% nz=(b_vector(:,3))*cos(theta)+(n_vector(:,3))*sin(theta);
% rx=bsxfun(@times,nx,radii);
% ry=bsxfun(@times,ny,radii);
% rz=bsxfun(@times,nz,radii*cline_para.zRatio);
% sx=bsxfun(@plus,CL2(:,1),rx);
% sy=bsxfun(@plus,CL2(:,2),ry);
% sz=bsxfun(@plus,size(im,3),rz);
% 
% 
%               maxR=sqrt(rkx.^2+rky.^2+rkz.^2);
%                 rkx=rkx./maxR;
%                 rkz=rkz./maxR;
%                 rky=rky./maxR;              
%               r_steps=linspace(min(maxR(:))*.8,1.2*max(maxR(:)),5);  %radius range
%                 
% %               r_steps=linspace(.8,1.2,5);
%                r_steps=reshape(r_steps,1,1,[]);  
%                 
%                 rkx=bsxfun(@times,rkx,r_steps);
%                 rky=bsxfun(@times,rky,r_steps);
%                 rkz=bsxfun(@times,rkz,r_steps);
%                 rkx=bsxfun(@plus,rx,CL3(:,1));
%                 rky=bsxfun(@plus,ry,CL3(:,2));
%                 rkz=bsxfun(@plus,rz,CL3(:,1)*0+midZ);
%                 
%                 
%                 
%                 
%                Ir=interp3(worm2,rkx,rky,rkz,'*linear');

% Is=squeeze((nanmax(Ir,[],2)));
% Is(isnan(Is))=0;
%%


%Is=smooth2a(Is,7,7);
%Is=imfilter(Is,gaussKernal2);
%imagesc((Is));
% V=permute(V,[3 2 1]);
% A=permute(A,[3 2 1]);
% 

V2=normalizeRange(V)-normalizeRange(A);
Vproj2=squeeze(nansum(V2,3));
Vproj=squeeze(nansum(V,3));

%%


wormProj=smooth2a(nanmax(V(:,:,:),[],3),1,1);
wormProj=normalizeRange(wormProj);
wormProjMask=wormProj>.2;

wormProj(wormProj<.2)=0;
%wormProj(:,1:end/2)=0;
wormProjMax=imregionalmax(wormProj);
wormProjMaxI=wormProj(wormProjMax);
%wormProjMaxI(wormProjMaxI>.8)=0;
[maxX,maxY]=find(wormProjMax);
[maxX,ia2]=sort(maxX);
maxY=maxY(ia2);
wormProjMaxI=wormProjMaxI(ia2);
maxPos=mean(maxX(wormProjMaxI>.8));
wormProjMaxI(maxX>(maxPos-20) & maxX<(maxPos+20))=0;
%%
[~,locs]=findpeaks(maxY);
locs=union(locs,find(maxY>max(maxY(locs))));
locs=union(locs,find(maxX>maxPos+20));
maxY=maxY(locs);
maxX=maxX(locs);
wormProjMaxI=wormProjMaxI(locs);


if 0% nnz(wormProjMaxI)>4;
f=fit(maxX,maxY,'poly2','Weights',max(0,(maxY-max(50,min(maxY))).*wormProjMaxI));
%[p,S,mu] = polyfit(maxX,maxY,2);
x=1:300;
yshift=f(x);
else
    yshift=ones(300,1);
end



%%
if counter~=1
    xIm=xcorr2(Vproj,Vtemplate);
[xlag,ylag]=find(xIm==max(xIm(:)));
lags=[xlag,ylag]-size(Vtemplate);
%[ismatrix(Vtemplate),ismatrix(Vproj)]


yshift=yshift-yshift(150+lags(1));
yshift(yshift>100)=100;
yshift(yshift<-100)=-100;
lags(2)=0;
%Vproj=movepixels(Vproj,lags(1),lags(2),0);
%Vproj=movepixels_2d_double(Vproj,lags(1),lags(2),0);

for i=1:size(V,3);
    V(:,:,i)=movepixels_2d_double(V(:,:,i),0,yshift,1);
    A(:,:,i)=movepixels_2d_double(A(:,:,i),0,yshift,1);

    V(:,:,i)=movepixels_2d_double(V(:,:,i),lags(1),lags(2),1);
    A(:,:,i)=movepixels_2d_double(A(:,:,i),lags(1),lags(2),1);
    
    
end


else
    lags=[0 0];
    Vtemplate=Vproj;
end

%% Transfrom control points

Fpoints2=zeros(size(ctrlPoints,1),3);
ctrlPoints2=cell2mat(ctrlPoints(:,[1 2 4]));
ctrlPoints2(:,3)=ctrlPoints2(:,3)-min(hiResIdx)+1;
for i=1:size(ctrlPoints,1);
Fpointsz=interp1(hiResIdx,1:length(hiResIdx),ctrlPoints{i,4},'nearest*','extrap');
Fpointsz(Fpointsz>length(hiResIdx))=length(hiResIdx);
Fpointsz(Fpointsz<1)=1;

Fpointsz=(Fpointsz-midZ)*3+window;
    xtemp=xslice(:,:,Fpointsz);
    ytemp=yslice(:,:,Fpointsz);
    xtemp(isnan(xtemp))=Inf;
        ytemp(isnan(xtemp))=Inf;
if ~isempty([ctrlPoints{i,1:2}])
    k=dsearchn([ytemp(:),xtemp(:)],[ctrlPoints{i,[2 1]}]);
    [Fpoints2(i,2),Fpoints2(i,1),Fpoints2(i,3)]=ind2sub(size(xtemp),k);
Fpoints2(i,3)=Fpointsz;
Fpoints2(i,1)=Fpoints2(i,1)-lags(2);
Fpoints2(i,2)=Fpoints2(i,2)-lags(1);
else
    Fpoints2(i,:)=[nan nan nan];
end
 Fpoints2(i,:)

    
end




%% 
xdip=104;

close all
V1=sum(V(:,:,1:100),3);
V1f=bpass(V1,.1,[30 30]);
V2=nansum(V(:,:,100:end),3);
V2f=bpass(V2,.1,[30 30]);


vline1=squeeze(sum(nansum(V1f,3),2));
vline2=squeeze(sum(nansum(V2f,3),2));
vline1=smooth(vline1,10);
vline2=smooth(vline2,10);
V1f=imfilter(V1f,-gaussKernal50);
V2f=imfilter(V2f,-gaussKernal50);
[~,maxpeak1]=max(vline1+vline2);
[~,maxpeak2]=max(vline2+vline1);

[~,locs1,w1,p1]=findpeaks(-vline1(1:maxpeak1));
w1=w1.*p1;
[~,nring1]=max(w1);
nring1=locs1(nring1);

[~,locs2,w2,p2]=findpeaks(-vline2(1:maxpeak2));
w2=w2.*p2;
[~,nring2]=max(w2);
nring2=locs2(nring2);
lowVal=min(V2f(:));
V2f(1:nring2-10,:)=lowVal;
V1f(1:nring1-10,:)=lowVal;


% V1f=[V1f; V1f*0-1000];
% V2f=[V2f; V2f*0-1000];

[~,locs3,w3,p3]=findpeaks(smooth(-vline1-vline2,15));
%if problem finding nchord peak, use old one
if any(locs3(locs3>maxpeak1))
nchord=locs3(locs3>maxpeak1);
nchord=nchord(1);


end
V2f(nchord:end,:)=lowVal;
V1f(nchord:end,:)=lowVal;
if counter==1
subV=V(nchord:end,:,:);
subV=squeeze(sum(subV,1));
[VXgrid,VYgrid]=meshgrid(1:size(subV,1),1:size(subV,2));
VXgrid=VXgrid-mean(VXgrid(:));
VYgrid=VYgrid-mean(VYgrid(:));
%Vbit=atan2(VXgrid,VYgrid);


Vbit=double(VXgrid>VYgrid)+2.*(double(VYgrid>-VXgrid));
Vbit=Vbit+1;
vRegion=accumarray(Vbit(:),subV(:),[],@nanmean);
vRegion=Vbit==find(vRegion==max(vRegion));
end
%%

xVec1=101:size(V1,2);
% aOut1(1,:)=fminsearch(@(a) CLsearch(V1f,xVec1,a(1)*(xVec1-101).^2+a(2)*(xVec1-101)+a(3) ...
%     ,show),[0.010 0.010 nring1]);
xVec2=1:100;
% aOut2(1,:)=fminsearch(@(a) CLsearch(V1f,xVec2,a(1)*(xVec2-101).^2+a(2)*(xVec2-101)+a(3) ...
%     ,show),[0.0040 0.0040 nring1]);
% 
% aOut1(2,:)=fminsearch(@(a) CLsearch(V2f,xVec1,a(1)*(xVec1-101).^2+a(2)*(xVec1-101)+a(3) ...
%     ,show),[0.0040 0.0040 nring2]);
% aOut2(2,:)=fminsearch(@(a) CLsearch(V2f,xVec2,a(1)*(xVec2-101).^2+a(2)*(xVec2-101)+a(3) ...
%     ,show),[0.0040 0.0040 nring2]);
% 
% boundary1=[aOut2(1,1)*(xVec2-101).^2+aOut2(1,2)*(xVec2-101)+aOut2(1,3) ...
%      aOut1(1,1)*(xVec1-101).^2+aOut1(1,2)*(xVec1-101)+aOut1(1,3)] ;
% boundary2=[aOut2(2,1)*(xVec2-101).^2+aOut2(2,2)*(xVec2-101)+aOut2(1,3) ...
%      aOut1(2,1)*(xVec1-101).^2+aOut1(2,2)*(xVec1-101)+aOut1(2,3);] ;

 
 aOut1=fminsearch(@(a) CLsearch(V1f,[xVec2 xVec1] ,[(a(1))*(xVec2-101).^2+a(2)*(xVec2-101)+nring1 ,...
    (a(3))*(xVec1-101).^2+a(4)*(xVec1-101)+nring1] ...
,show),[0.0040 0.0040  .004 .3]);

%  aOut1=fminsearch(@(a) CLsearch(V1f,[xVec2 xVec1] ,[a(4)*(sinh(a(1)*(xVec2-101)))+a(2) ,...
%     a(5)*(sinh(a(3)*(xVec1-101)))+a(2)] ...
% ,show),[0.040 nring1 .04 1 1]);
% 

 aOut2=fminsearch(@(a) CLsearch(V2f,[xVec2 xVec1] ,[(a(1))*(xVec2-101).^2+a(2)*(xVec2-101)+nring2 ,...
   ( a(3))*(xVec1-101).^2+a(4)*(xVec1-101)+nring2] ...
,show),[0.0040 -0.0040  .004 .3  ]);


boundary1=[aOut2(1)*(xVec2-101).^2+aOut2(2)*(xVec2-101)+nring1 ...
     aOut1(3)*(xVec1-101).^2+aOut1(4)*(xVec1-101)+nring1] ;
boundary2=[aOut2(1)*(xVec2-101).^2+aOut2(2)*(xVec2-101)+nring2 ...
     aOut2(3)*(xVec1-101).^2+aOut2(4)*(xVec1-101)+nring2] ;

 
 boundary2(boundary2>nchord)=nchord;
  boundary1(boundary1>nchord)=nchord;


 
 %% cut up regions
 wormRegions=zeros(size(V));
 [wormy,wormx]=meshgrid(1:size(V,2),1:size(V,1));
wormRegions1=bsxfun(@ge,wormx,boundary1); 
wormRegions2=bsxfun(@ge,wormx,boundary2); 
wormRegions(:,:,1:100)=repmat(wormRegions1,1,1,100);
wormRegions(:,:,101:end)=repmat(wormRegions2,1,1,101);

wormRegions(nchord:end,:,:)=wormRegions(nchord:end,:,:)+...
2*permute(repmat(vRegion, 1,1,1+size(V,1)-nchord),[3,1,2]);
wormRegions(wormRegions==2)=0;
wormRegions(wormRegions==3)=2;


%%
%CLsearch(V1,xVec,a*(xVec-101).^2+104,show)

% Vall{counter}=V;
% Aall{counter}=A;
% IsAll{counter}=Is;
% IrAll{counter}=Ir;
% rawAll{counter}=worm2;
% VprojAll{counter}=Vproj;
% VprojAll2{counter}=Vproj;
fileName=[imageFolder filesep 'image' num2str(iStack,'%3.5d') '.tif'];
fileName2=[imageFolder2 filesep 'image' num2str(iStack,'%3.5d') '.tif'];
fileName3=[imageFolder2 filesep 'imageMap' num2str(iStack,'%3.5d') '.tif'];

%tiffwrite(fileName,Vproj,'tif');
tiffwrite(fileName2,V,'tif');
tiffwrite(fileName3,wormRegions,'tif');
fpointsAll{iStack}=Fpoints2;
display(['Finished image ' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])
iStack
fclose(Fid);

     catch ME
        ME
        display(['Error in Frame' num2str(iStack,'%3.5d') ' in ' num2str(toc) 's'])

    end
    
    
end
%%
for iStack=startStack:endStack
    fileName4=[imageFolder2 filesep 'controlPoints' num2str(iStack,'%3.5d')];
    Fpoints2=fpointsAll{iStack};
    save(fileName4,'Fpoints2');
end






%%

for iStack=startStack:endStack
fileName2=[imageFolder2 filesep 'image' num2str(iStack,'%3.5d') '.tif'];
fileName3=[imageFolder2 filesep 'imageMap' num2str(iStack,'%3.5d') '.tif'];
wormRegions=stackLoad(fileName3);
worm=stackLoad(fileName2);
    
    
end




%%
[Fy,Fx,Fz]=gradient(im);




fx=interp3(Fx,sx,sy,sz,'*linear');
fy=interp3(Fy,sx,sy,sz,'*linear');
fz=interp3(Fz,sx,sy,sz,'*linear');

fr=nx.*fx+ny.*fy+nz.*fz;
fr=nansum(fr,2);
fx=nansum(fx,2);
fy=nansum(fy,2);
fz=nansum(fz,2);

%%








%wormSmooth=bpass3(worm,10,[20,20,15]);



