%This is Jeff's amazing code!! We are so greatfule for jeff 
% To work on it

% %% load George's file with manual identification
% A=importdata('/Users/jeffnguyen/Documents/Data/for andy/d1.tif.marker');
% A.x=A.data(:,1);
% A.y=A.data(:,2);
% A.z=A.data(:,3)*2;
% 

%% load worm image


worm=stackLoad('C:\Users\gplummer\Documents\data\odr2\03242014\odr2_gfpa2.tif', 34,2);
imsize=size(worm);

%% Initialize parameters
thresh1=.03; %initial Threshold
hthresh=-.0001; %threshold for trace of hessian.
minObjSize=100; 
maxObjSize=1500;
minObjectSpacing=5;
minSearchRad=3;
pad=4;
show=1;
box=true(3,3,3);

%% subtract pedistal, normalize, filter
pedMask=false(imsize);
pedMask(1:3,:,:)=true;
pedMask(:,1:3,:)=true;
pedMask(end-2:end,:,:)=true;
pedMask(:,end-2:end,:)=true;

pedistal=median(worm(pedMask));
worm=worm-pedistal;
wormtop=imtophat(worm,strel('ball',25,25,0)); %top hat filter (SLOW!!!)
wormtop=worm-abs(worm-wormtop);
wormtop(wormtop<0)=0;
wormtop=normalizeRange(wormtop);

%% initial threshold
wormBW=wormtop>thresh1;
%% find connected objects
cc=bwconncomp(wormBW,6);
blobSizes=cellfun(@(x) length(x), cc.PixelIdxList);
cc.PixelIdxList(blobSizes<minObjSize)=[];
cc.NumObjects=sum(blobSizes>minObjSize);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



%% use hessian to find nuclei in objects
centerIm=zeros(size(wormBW));
wormBW2=zeros(size(wormBW));
for iblob=1:cc.NumObjects;
    %crop out object with pad
    BB=floor(blobStats(iblob).BoundingBox);
    BB(1:3)=BB(1:3)-[pad,pad,pad];
    BB(4:end)=BB(4:end)+BB(1:3)+[2*pad,2*pad,2*pad];
    %don't overshoot size of image
    overEdge=[1 1 1 imsize]-BB;
    overEdge(4:6)=-overEdge(4:6);
    overEdge=max(overEdge,zeros(1,6));   
    BB(BB<1)=1;
    BB([false,false,false,bsxfun(@ge,BB(4:6),imsize)])=imsize((bsxfun(@ge,BB(4:6),imsize)));


    subBW=wormBW(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    subIm=wormtop(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    
    %filter/normalize sub image
    if mean(size(subIm))<30
        subIm=imtophat(subIm,strel('disk',10'));
    end
    
    subIm=normalizeRange(subIm);
%     minIm=imerode(subIm,true(11,11,1));
%     maxIm=imdilate(subIm,true(11,11,1));
%    maxIm=smooth3(maxIm,'box',[5,5,5]);
%     
%     subThresh=(subIm-minIm);
%     subThresh(maxIm<.01)=0;
%     subThresh=subThresh./(maxIm-minIm);
%     
%subBW=subIm>.1;

pedMask=true(size(subIm));
pedMask(1:pad-overEdge(2),:,:)=false;
pedMask(:,1:pad-overEdge(1),:)=false;
pedMask(:,:,1:pad-overEdge(3))=false;
pedMask(end-pad+1+overEdge(5):end,:,:)=false;
pedMask(:,end-pad+1+overEdge(4):end,:)=false;
pedMask(:,:,end-pad+1+overEdge(6):end)=false;
  
    

% BWe=bwulterode(subBW);
% [x,y,z]=ind2sub(size(subBW),find(BWe));
% %pixel watershed
% DD=bwdist(~subBW);
% DD=-DD;
% DD(~subBW)=Inf;
% DD=imhmin(DD,.5);
% DL=watershed(DD,18);

% smooth image and calculate hessian and eigenvalues
subIm=smooth3(subIm,'gaussian',5,3);
H=hessianMatrix(subIm,3);
Heig=hessianEig(H,subBW);
Htrace=sum(Heig,4);
% Jm= Heig(:,:,:,1)<-hthresh & Heig(:,:,:,2)<-hthresh & ...
%    Heig(:,:,:,3)<-hthresh;
Jm=Htrace<hthresh;
Jm=Jm & pedMask;
Jm=xyzConvHull(Jm,3);

% watershed filter shapes
Jd=-bwdist(~Jm);
%Jd=smooth3(Jd,'gaussian',5,2);
Jd=imhmin(Jd,.8);
Jd(~Jm)=Inf;
Jw=watershed(Jd);
Jm=Jm.*(Jw>0);
Jm=imerode(Jm,true(2,2,2));

subCC=bwconncomp(Jm>0,6);

%check size, if size is too large, increase watershedding
for isubBlob=1:subCC.NumObjects
    if length(subCC.PixelIdxList{isubBlob})>50
        blank=false(size(subBW));
        blank(subCC.PixelIdxList{isubBlob})=true;
        [x,y,z]=ind2sub(size(subBW),subCC.PixelIdxList{isubBlob});
 %       [coeff,score,latent,~,explianed]=pca([x,y,z]);
%        sp=sum((latent-circshift(latent,1)).^2)./sum(latent.^2);
      %  Jm(subCC.PixelIdxList{isubBlob})=length(x);
        if length(x)>maxObjSize
            Jw=ones(size(Jd));
            watershedthresh=.7;
            while(all(Jw(:))) || ~sum(~Jw(blank))
            Jd=-bwdist(~blank);
            %Jd=smooth3(Jd,'gaussian',5,2);
            Jd=imhmin(Jd,watershedthresh);
            Jd(~Jm)=Inf;
            Jw=watershed(Jd);
            watershedthresh=watershedthresh-.1;
            end
            Jm(blank)=~~Jw(blank);
            
        end
    end
    
end

% remove small objects in subImage
subCC=bwconncomp(Jm>0,6);
subBlobStats=regionprops(subCC,'Centroid');
objectSizes=cellfun(@(x) length(x), subCC.PixelIdxList);
badIdx=objectSizes<minObjSize;
bad=cell2mat(subCC.PixelIdxList(badIdx)');
subCC.PixelIdxList((badIdx))=[];
subCC.NumObjects=sum(~badIdx);
Jm(bad)=0;

centerPnts=round([subBlobStats.Centroid]);
centerPnts=reshape(centerPnts',3,[])';
Jc=false(size(Jm));
Jc(sub2ind(size(Jm),centerPnts(:,2),centerPnts(:,1),centerPnts(:,3)))=true;
Jc(bad)=false;




%% find centroids, display
if show
centerIm(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6))=Jc;
imagesc(sum(subIm,3));
hold on
axis equal
temp=sum(Jc,3);
[y,x]=find(temp);
scatter(x,y,'black.');
hold off
pause(1);
end

%compile results of all sub images into final segmented mask. 
wormBW2(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6))=Jm;
% 
% 
%         
        
end
%% 
%centerIm is binary image of all centroid positions. 
centerIm=bwulterode(centerIm);

[y,x,z]=ind2sub(size(wormBW),find(centerIm));
    figure;
scatter3(x,y,z);axis equal
% %% show George's result as well. 
% hold on
% 
% scatter3(A.x,A.y,A.z+mean(z)-mean(A.z),'rx');
% 
% disp([ num2str(length(x)) ' points found, compared to ' num2str(length(A.x)) ,...
%     ' found by Geroge']);
% 
% 

