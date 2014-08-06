function wormBW2=WormSegmentHessian3dwhole(worm,options)

%% Initialize parameters
thresh1=.03; %initial Threshold
hthresh=-.0001; %threshold for trace of hessian.
minObjSize=500; 
maxObjSize=Inf;
watershedFilter=0;
filterSize=[50,50,50];
pad=4;
noise=9;
show=0;
if nargin==2
    Fnames=fieldnames(options);
    for i=1:length(Fnames)
        eval([Fnames{i} '= options.' Fnames{i} ';']);
    end
    
end




imsize=size(worm);
imsize=imsize([2,1,3]);
%% subtract pedistal, normalize, filter

pedMask=false(imsize);
pedMask(1:3,:,:)=true;
pedMask(:,1:3,:)=true;
pedMask(end-2:end,:,:)=true;
pedMask(:,end-2:end,:)=true;

pedistal=median(worm(pedMask));
worm=worm-pedistal;
worm(worm<0)=0;

%wormtop=imtophat(worm,strel('ball',25,25,0)); %top hat filter (SLOW!!!)
wormtop=bpass3_jn(worm,noise,filterSize);
wormtop=smooth3(worm,'gaussian',[11,11,11],noise);

%% initial threshold
wormBW=wormtop>thresh1;
%% find connected objects
% cc=bwconncomp(wormBW,6);
% blobSizes=cellfun(@(x) length(x), cc.PixelIdxList);
% cc.PixelIdxList(blobSizes<minObjSize)=[];
% cc.NumObjects=sum(blobSizes>=minObjSize);
% blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



% smooth image and calculate hessian and eigenvalues
Jm=eigenSegmentation3d(wormtop,0);
Jm=xyzConvHull(Jm,1:3);
Jm=xyzConvHull(Jm,1:3);

% watershed filter shapes
if watershedFilter
Jd=-bwdist(~Jm);
%Jd=smooth3(Jd,'gaussian',5,2);
Jd=imhmin(Jd,watershedFilter);
Jd(~Jm)=Inf;
Jw=watershed(Jd);
Jm=Jm.*(Jw>0);
end
Jm=imerode(Jm,true(2,2,2));



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
            while((all(Jw(:))) || ~sum(~Jw(blank))) && watershedthresh>.4
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
        
%% 
%centerIm is binary image of all centroid positions. 
centerIm=bwulterode(centerIm);

%[y,x,z]=ind2sub(size(wormBW),find(centerIm));
%    figure;
%scatter3(x,y,z);axis equal

