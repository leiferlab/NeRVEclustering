function wormBW2=WormSegmentHessian2D_whole(worm)

%% Initialize parameters
thresh1=.05; %initial Threshold
hthresh=-.0001; %threshold for trace of hessian.
minObjSize=20; 
maxObjSize=200;
minObjectSpacing=5;
minSearchRad=3;
pad=4;
show=0;
watershedthresh=0;
maxSplit=0;

box=true(3,3);
imsize=size(worm);
imsize=imsize([2,1]);
%% subtract pedistal, normalize, filter

pedMask=false(imsize);
pedMask(1:3,:)=true;
pedMask(:,1:3)=true;
pedMask(end-2:end,:)=true;
pedMask(:,end-2:end)=true;

pedistal=nanmedian(worm(pedMask));
worm=worm-pedistal;
worm(worm<0)=0;
worm(isnan(worm))=0;
%wormtop=imtophat(worm,strel('ball',25,25,0)); %top hat filter (SLOW!!!)
wormtop=bpass_jn(worm,3,[30,30]);
%wormtop=worm-abs(worm-wormtop);
%wormtop(wormtop<0)=0;
wormtop=normalizeRange(wormtop);

%% initial threshold
wormBW=wormtop>thresh1;
%% find connected objects
cc=bwconncomp(wormBW,4);
blobSizes=cellfun(@(x) length(x), cc.PixelIdxList); % replace with bwareaopen
cc.PixelIdxList(blobSizes<minObjSize)=[];
cc.NumObjects=sum(blobSizes>=minObjSize);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



%% use hessian to find nuclei in objects
centerIm=zeros(size(wormBW));
wormBW2=zeros(size(wormBW));

% smooth image and calculate hessian and eigenvalues
% subIm=normalizeRange(smooth2a(wormtop,10,10));
% subIm=subIm./smooth2a(subIm,30,30).*wormBW;
% 
subIm=wormtop;
H=hessianMatrix(subIm,5);
Heig=hessianEig(H,wormBW);

Htrace=sum(Heig,3);

Jm=Htrace<hthresh;
%Jm=((Jm | subIm>.2) & wormtop>thresh1); %add global thresh
%%
Jm=AreaFilter(Jm,minObjSize,Inf,4);
subCC=bwconncomp(Jm>0,4);
stats=regionprops(subCC,'eccentricity','Area');
%%

%check size, if size is too large, increase watershedding
% for isubBlob=1:subCC.NumObjects
%     if 1% stats(isubBlob).Eccentricity>.8 && stats(isubBlob).Area>maxObjSize;
%         blank=false(size(worm));
%         blank(subCC.PixelIdxList{isubBlob})=true;
%             Jd=-bwdist(~blank);
%             %Jd=smooth3(Jd,'gaussian',5,2);
%             Jd=imhmin(Jd,watershedthresh);
%             Jd(~Jm)=Inf;
%             Jw=watershed(Jd);
%           %  Jw=watershed(-imregionalmax(blank.*wormtop));
%             
%             
%             Jm(blank)=~~Jw(blank);
%     end
%     
% end

   Jd=-bwdist(~Jm);
            %Jd=smooth3(Jd,'gaussian',5,2);
            Jd=imhmin(Jd,watershedthresh);
            Jd(~Jm)=Inf;
            Jw=watershed(Jd);
          %  Jw=watershed(-imregionalmax(blank.*wormtop));
            
            Jm=Jm.*(Jw>0);
            %Jm(blank)=~~Jw(blank);

%%
if maxSplit
    subCC=bwconncomp(Jm>0,4);
%stats=regionprops(subCC,'eccentricity','Area');


%check size, if size is too large, increase watershedding
for isubBlob=1:subCC.NumObjects
  %  if stats(isubBlob).Eccentricity>.8 && stats(isubBlob).Area>maxObjSize;
        blank=false(size(worm));
        blank(subCC.PixelIdxList{isubBlob})=true;
            Jw=watershed(-imregionalmax(blank.*wormtop));
            
            
            Jm(blank)=~~Jw(blank);
%    end
    
end
    
    
end


% remove small objects in subImage
% subCC=bwconncomp(Jm>0,4);
% subBlobStats=regionprops(subCC,'Centroid');
% objectSizes=cellfun(@(x) length(x), subCC.PixelIdxList);
% badIdx=objectSizes<minObjSize;
% bad=cell2mat(subCC.PixelIdxList(badIdx)');
% subCC.PixelIdxList((badIdx))=[];
% subCC.NumObjects=sum(~badIdx);
% Jm(bad)=0;
% 
% centerPnts=round([subBlobStats.Centroid]);
% centerPnts=reshape(centerPnts',2,[])';
% Jc=false(size(Jm));
% Jc(sub2ind(size(Jm),centerPnts(:,2),centerPnts(:,1)))=true;
% Jc(bad)=false;
%
%compile results of all sub images into final segmented mask.
Jm=AreaFilter(Jm,minObjSize,Inf,4);

wormBW2=Jm;
% 
% 
%         
        
%% 
%centerIm is binary image of all centroid positions. 
%centerIm=bwulterode(centerIm);

%[y,x,z]=ind2sub(size(wormBW),find(centerIm));
%    figure;
%scatter3(x,y,z);axis equal

