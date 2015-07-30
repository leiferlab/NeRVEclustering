function [wormBW2,wormtop]=WormSegmentHessian3dStraighten(worm,options)

%% Initialize default parameters, all of these can also be fields in options
thresh1=.03; %initial Threshold
hthresh=-.0; %threshold for trace of hessian.
minObjSize=80; % min object size
maxObjSize=5000; % max object size
valleyRatio=.75;
watershedFilter=0; % watershed filter object shapes? is also the value for imhmin
filterSize=[10,10,4]; %bp filter size low f
noise=1; % bp filter hi f
pad=10; % pad to take around each sub blob
show=0; %show fits (deactivated)
maxSplit=0; % split objects using regional maxima
minSphericity=.55; % minimum sphericity for splitting.
prefilter=0;
gaussFilter=1;

% parse options to load fields
if nargin==2
    Fnames=fieldnames(options);
    for i=1:length(Fnames)
        eval([Fnames{i} '= options.' Fnames{i} ';']);
    end
    
end

%cuberoot for min obj dimension
minObjDim=round(minObjSize^(1/3));
imsize=size(worm);
imsize=imsize([2,1,3]);
%% subtract pedistal, normalize, filter

worm=pedistalSubtract(worm);
worm(worm<0)=0;

if ~prefilter
%wormtop=imtophat(worm,strel('ball',25,25,0)); %top hat filter (SLOW!!!)
wormtop=bpass3(worm,noise,filterSize);
else
    wormtop=worm;
end

wormtop=normalizeRange(wormtop);

%% initial threshold
wormBW=wormtop>thresh1;
%wormBW=imclearborder(wormBW,6);
wormBW=AreaFilter(wormBW,minObjSize,[],6);

%% find connected objects
cc=bwconncomp(wormBW,6);
blobSizes=cellfun(@(x) length(x), cc.PixelIdxList);
cc.NumObjects=sum(blobSizes>=minObjSize);
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
%     if mean(size(subIm))<30
%         subIm=imtophat(subIm,strel('disk',10'));
%     end
    
    subIm=normalizeRange(subIm);


% pedMask=true(size(subIm));
% pedMask(1:pad-overEdge(2),:,:)=false;
% pedMask(:,1:pad-overEdge(1),:)=false;
% pedMask(:,:,1:pad-overEdge(3))=false;
% pedMask(end-pad+1+overEdge(5):end,:,:)=false;
% pedMask(:,end-pad+1+overEdge(4):end,:)=false;
% pedMask(:,:,end-pad+1+overEdge(6):end)=false;
%   
%     

% BWe=bwulterode(subBW);
% [x,y,z]=ind2sub(size(subBW),find(BWe));
% %pixel watershed
% DD=bwdist(~subBW);
% DD=-DD;
% DD(~subBW)=Inf;
% DD=imhmin(DD,.5);
% DL=watershed(DD,18);

% smooth image and calculate hessian and eigenvalues for segmentation

 if ~prefilter
 subIm=smooth3(subIm,'gaussian',2*gaussFilter+1,gaussFilter);
 else
 subIm=bpass3(subIm,2,filterSize);
 end
H=hessianMatrix(subIm,8);
Heig=hessianEig(H,subBW);
Heig(isnan(Heig))=0;
Htrace=real(Heig(:,:,:,1));
% Jm= Heig(:,:,:,1)<-hthresh & Heig(:,:,:,2)<-hthresh & ...
%    Heig(:,:,:,3)<-hthresh;
Jm=Htrace<hthresh ;
Jm=AreaFilter(Jm,minObjSize,[],6);

%Jm=Jm & pedMask;
% watershed filter shapes
%%
if watershedFilter
Jd=-bwdist(~Jm);  %make distance map
%Jd=smooth3(Jd,'gaussian',5,2);
Jd=imhmin(Jd,watershedFilter);
Jd(~Jm)=Inf;
Jw=watershed(Jd);
Jm=Jm.*(Jw>0);
end

%%
%Jm=imclearborder(Jm);


%watershed splitting based on local maxima locations
if maxSplit
    %find regionalmaxima and threshold around that intensity
    subImaxPnts=imregionalmax(subIm.*Jm);


    %subImaxReg=subIm>(max(subIm(subImaxPnts)))*valleyRatio;
    subImax=imdilate(subImaxPnts.*Jm.*subIm*valleyRatio,true(minObjDim,minObjDim,minObjDim));
   subImaxReg=subImaxPnts>subImax & subIm>(2*thresh1);
    % subImax=imregionalmax(subIm);
    % 
    % subImax=imdilate(subImax.*Jm,true(minObjDim,minObjDim,minObjDim));

    JmLabel=bwlabeln(Jm,6);
    for iLabel=1:max(JmLabel(:));
        subJm=JmLabel==iLabel;
        subsubImax=subImaxReg & subJm;
        subsubcc=bwconncomp(subsubImax);
        if subsubcc.NumObjects>1
            maxBW=watershed(bwdist(subsubImax));

            Jm(maxBW==0 & subJm)=0;
        end

    end

end
%%
Jm=Jm.*subBW;
Jm=AreaFilter(Jm,minObjSize,[],6);

subCC=bwconncomp(Jm>0,6);

%check size, if size is too large, increase watershedding
for isubBlob=1:subCC.NumObjects
    if length(subCC.PixelIdxList{isubBlob})>minObjSize
        %calculate volume and sphericity for checking before further
        %watershed
        blank=false(size(subBW));
        blank(subCC.PixelIdxList{isubBlob})=true;
        [blanky,blankx,blankz]=ind2sub(size(blank),subCC.PixelIdxList{isubBlob});
        
 %       [x,y,z]=ind2sub(size(subBW),subCC.PixelIdxList{isubBlob});
        subblobP=bwperim(blank(min(blanky):max(blanky),...
            min(blankx):max(blankx),min(blankz):max(blankz)),26);
        subblobP=sum(subblobP(:));
        
        subblobV=sum(blank(:));
        subblobS=(pi*36*subblobV.^2)^(1/3)/subblobP;
        
        
        
 %       [coeff,score,latent,~,explianed]=pca([x,y,z]);
%        sp=sum((latent-circshift(latent,1)).^2)./sum(latent.^2);
      %  Jm(subCC.PixelIdxList{isubBlob})=length(x);
        if subblobV>maxObjSize || subblobS<minSphericity
%            [subblobV>maxObjSize , subblobS<minSphericity]
            Jw=ones(size(Jm));
            watershedthresh=.7;
            while((all(Jw(:))) || ~sum(~Jw(blank))) && watershedthresh>.4
            Jd=-bwdist(~blank);
            %Jd=smooth3(Jd,'gaussian',5,2);
            Jd=imhmin(Jd,watershedthresh);
            Jd(~blank)=Inf;
            subsubcc=bwconncomp(imregionalmin(Jd));
            if subsubcc.NumObjects>1
            Jw=watershed(Jd);
            end

            watershedthresh=watershedthresh-.1;
            end
            Jm(blank)=~~Jw(blank);
            
        end
    end
    
end

% remove small objects in subImage
Jm=Jm.*subBW;
Jm=AreaFilter(Jm,minObjSize,[],6);
%Jm=xyzConvHull(Jm,3);
% centerPnts=round([subBlobStats.Centroid]);
% centerPnts=reshape(centerPnts',3,[])';
% Jc=false(size(Jm));
% Jc(sub2ind(size(Jm),centerPnts(:,2),centerPnts(:,1),centerPnts(:,3)))=true;
% Jc(bad)=false;




%% find centroids, display (off)
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
%Jm=imclearborder(Jm);

wormBW2(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6))=or(Jm,wormBW2(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6)));
% 
% 
%         
        
end
%% 


%centerIm is binary image of all centroid positions. 

%[y,x,z]=ind2sub(size(wormBW),find(centerIm));
%    figure;
%scatter3(x,y,z);axis equal

