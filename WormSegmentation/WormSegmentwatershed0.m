
%% load worm image
worm=stackLoad('/Users/jeffnguyen/Documents/Data/for andy/d2.tif', 28,2);
A=importdata('/Users/jeffnguyen/Documents/Data/for andy/d2.tif.marker');
A.x=A.data(:,1);
A.y=A.data(:,2);
A.z=A.data(:,3)*2;

imsize=size(worm);

%%
thresh1=.03;
thresh2=.1;
minObjSize=100;
stdFiltRadius=5;
minObjectSpacing=5;
minSearchRad=3;
pad=4;

s2=stdFiltRadius*2-1;
ball=zeros(s2,s2,s2);
ball(stdFiltRadius,stdFiltRadius,stdFiltRadius)=1;
ball=bwdist(ball);
ball=ball<stdFiltRadius;
ball=smooth3(ball,'box',[3,3,3]);
box=true(3,3,3);

msk=true(minSearchRad,minSearchRad,minSearchRad);
msk(minSearchRad/2+.5,minSearchRad/2+.5,minSearchRad/2+.5)=false;


% msk=true(minSearchRad,minSearchRad,1);
% msk(minSearchRad/2+.5,minSearchRad/2+.5,:)=false;
% msk1=permute(msk,[2,3,1]);
% msk2=permute(msk,[3,1,2]);

%% subtract pedistal
pedMask=false(imsize);
pedMask(1:3,:,:)=true;
pedMask(:,1:3,:)=true;
pedMask(end-2:end,:,:)=true;
pedMask(:,end-2:end,:)=true;

pedistal=median(worm(pedMask));
worm=worm-pedistal;
wormtop=imtophat(worm,strel('ball',25,25,0));
wormtop=worm-abs(worm-wormtop);
wormtop(wormtop<0)=0;
wormtop=normalizeRange(wormtop);

%%
wormBW=wormtop>thresh1;
%%
cc=bwconncomp(wormBW,6);
blobSizes=cellfun(@(x) length(x), cc.PixelIdxList);
cc.PixelIdxList(blobSizes<minObjSize)=[];
cc.NumObjects=sum(blobSizes>minObjSize);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



%%
centerIm=zeros(size(wormBW));

for iblob=1:cc.NumObjects;
    BB=floor(blobStats(iblob).BoundingBox);
    BB(1:3)=BB(1:3)-[pad,pad,pad];
    BB(4:end)=BB(4:end)+BB(1:3)+[2*pad,2*pad,2*pad];
    overEdge=[1 1 1 imsize]-BB;
    overEdge(4:6)=-overEdge(4:6);
    overEdge=max(overEdge,zeros(1,6));
    
    BB(BB<1)=1;
    BB([false,false,false,bsxfun(@ge,BB(4:6),imsize)])=imsize((bsxfun(@ge,BB(4:6),imsize)));


    
    subBW=wormBW(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    subIm=wormtop(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    if mean(size(subIm))<30
        subIm=imtophat(subIm,strel('disk',10'));
    end
    
    subIm=normalizeRange(subIm);
    minIm=imerode(subIm,true(11,11,1));
    maxIm=imdilate(subIm,true(11,11,1));
   maxIm=smooth3(maxIm,'box',[5,5,5]);
    
    subThresh=(subIm-minIm);
    subThresh(maxIm<.01)=0;
    subThresh=subThresh./(maxIm-minIm);
    
%subBW=subIm>.1;
    pedMask=true(size(subIm));
pedMask(1:pad-overEdge(1),:,:)=false;
pedMask(:,1:pad-overEdge(2),:)=false;
pedMask(:,:,1:pad-overEdge(3))=false;
pedMask(end-pad+1+overEdge(4):end,:,:)=false;
pedMask(:,end-pad+1+overEdge(5):end,:)=false;
pedMask(:,:,end-pad+1+overEdge(6):end)=false;
  
    

% BWe=bwulterode(subBW);
% [x,y,z]=ind2sub(size(subBW),find(BWe));
% %pixel watershed
% DD=bwdist(~subBW);
% DD=-DD;
% DD(~subBW)=Inf;
% DD=imhmin(DD,.5);
% DL=watershed(DD,18);

% std filt methodc
J=smooth3(subIm,'gaussian',5,3);
J=stdfiltN(J,ball);
J=normalizeRange(J);
J=imhmin(normalizeRange(J),.0001);
J=-J;
Jm=J>imdilate(J,msk);
%Jm=imdilate(J>Jm1,box)+imdilate(J>Jm2,box)+imdilate(J>Jm3,box);
%Jm=(J>Jm1)+(J>Jm2)+(J>Jm3);
%Jm=Jm>1;
Jm=Jm.*subBW;
%Jm=imclearborder(Jm);


Jm=imdilate(Jm,true(minObjectSpacing,minObjectSpacing,minObjectSpacing));
Jm=Jm.*(subThresh>.4);
Jm=Jm.*pedMask;
subCC=bwconncomp(Jm>0,6);
subBlobStats=regionprops(subCC,'Centroid');
centerPnts=round([subBlobStats.Centroid]);
centerPnts=reshape(centerPnts',3,[])';
Jm=zeros(size(Jm));
Jm(sub2ind(size(Jm),centerPnts(:,2),centerPnts(:,1),centerPnts(:,3)))=1;

sum(Jm(:))
centerIm(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6))=Jm;
imagesc(sum(subIm,3));
hold on
axis equal
temp=sum(Jm,3);
[y,x]=find(temp);
scatter(x,y,'black.');
hold off
pause(1);
% 
% subBW2=subBW.*(DL>0 );%| subIm>thresh2);
% subIm2=subIm.*(DL>0 | subIm>thresh2).*subBW;
%     subCC=bwconncomp(subBW2,18);
%     objectSizes=cellfun(@(x) length(x), subCC.PixelIdxList);
% bad=cell2mat(subCC.PixelIdxList(objectSizes<100)');
% subBW2(bad)=0;
% wormBW(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6))=subBW2;
% 
% 
%         
        
end
%%
centerIm=bwulterode(centerIm);
[y,x,z]=ind2sub(size(wormBW),find(centerIm));
    figure;
scatter3(x,y,z);axis equal
hold on

scatter3(A.x,A.y,A.z+mean(z)-mean(A.z),'rx');

disp([ num2str(length(x)) ' points found, compared to ' num2str(length(A.x)) ,...
    ' found by Geroge']);
