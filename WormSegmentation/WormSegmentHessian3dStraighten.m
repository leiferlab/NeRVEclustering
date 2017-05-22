function [im_bw_out,im_smooth]=...
    WormSegmentHessian3dStraighten(im,options,im_smooth)
% this function takes as an input an image with bright blobs and segments
% it by looking at the eigenvalues of the hessian matrix. Objects can be
% further seperated using a watershed filter based on the object size. An
% optional 3rd input is a smoothed version of the image, which bypasses the
% need to do a bpass filter in this code. 

%% Initialize default parameters, all of these can also be fields in options
thresh1=.03; %initial Threshold
hthresh=-.0; %threshold for trace of hessian.
minObjSize=80; % min object size
maxObjSize=350; % max object size
valleyRatio=.75;
% watershed filter object shapes? is also the value for imhmin
watershedFilter=0; 
filterSize=[10,10,4]; %bp filter size low f
noise=1; % bp filter hi f
pad=10; % pad to take around each sub blob
show=0; %show fits (deactivated)
maxSplit=0; % split objects using regional maxima
minSphericity=.84; % minimum sphericity for splitting.
prefilter=0;
gaussFilter=1;

% parse options to load fields
if nargin>=2
    Fnames=fieldnames(options);
    for i=1:length(Fnames)
        eval([Fnames{i} '= options.' Fnames{i} ';']);
    end
    
end

%cuberoot for min obj dimension
minObjDim=round(minObjSize^(1/3));
imsize=size(im);
imsize=imsize([2,1,3]);
%% if smoothed image is an input, use it, otherwise, do a bpass filter

im(im<0)=0;
if ~prefilter
im_smooth=bpass3(im,noise,filterSize);
else
    if nargin<3
    im_smooth=im;
    end
end

im_smooth=normalizeRange(im_smooth);

%% initial threshold
im_bw=im_smooth>max(graythresh(im_smooth(:)),thresh1);
%remove small objects
im_bw=AreaFilter(im_bw,minObjSize,[],6);

% find connected objects
cc=bwconncomp(im_bw,6);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');



%% use hessian to find nuclei in objects
% after initial rough thresholding, loop through each segmented region and
% do a thresholding based on hessian and watersheds. This is faster than
% just shotgun analyzing the entire image.


im_bw_out=zeros(size(im_bw));
for iblob=1:cc.NumObjects
    %crop out object with pad
    box=floor(blobStats(iblob).BoundingBox);
    box(1:3)=box(1:3)-[pad,pad,pad];
    box(4:end)=box(4:end)+box(1:3)+[2*pad,2*pad,2*pad];
    
    %don't overshoot size of image if negative, cast to 0, if larger than
    %imsize, make the edge of the box the edge of the image. 
    box(box<1)=1;
    over_edge=bsxfun(@ge,box(4:6),imsize);
    box([false,false,false,over_edge])=imsize(over_edge);


    sub_bw=im_bw(box(2):box(5),box(1):box(4),box(3):box(6));
    sub_im=im_smooth(box(2):box(5),box(1):box(4),box(3):box(6));
    sub_im=normalizeRange(sub_im);

% smooth image and calculate hessian and eigenvalues for segmentation,
% filter less if prefiltered
 if ~prefilter
 sub_im=smooth3(sub_im,'gaussian',2*gaussFilter+1,gaussFilter);
 else
 sub_im=bpass3(sub_im,2,filterSize);
 end
%clculate hessian matrix for each point
H=hessianMatrix(sub_im,8);
%Find the eigenvalues for the hessian at each point which is above
%threshold. 
Heig=hessianEig(H,sub_bw);
Heig(isnan(Heig))=0;
%find where the trace is above a threshold
Htrace=real(Heig(:,:,:,1));
hess_bw=Htrace<hthresh ;
%apply area threshold
hess_bw=AreaFilter(hess_bw,minObjSize,[],6);

%% watershed filter shapes

if watershedFilter
Jd=-bwdist(~hess_bw);  %make distance map
%Jd=smooth3(Jd,'gaussian',5,2);
Jd=imhmin(Jd,watershedFilter);
Jd(~hess_bw)=Inf;
Jw=watershed(Jd);
hess_bw=hess_bw.*(Jw>0);
end

%%
%watershed splitting based on local maxima locations
if maxSplit
    %find regionalmaxima and threshold around that intensity
    subImaxPnts=imregionalmax(sub_im.*hess_bw);
    
    %make regions around each maxima that was segmented by the hessian.
    %The regions have value of the intensity of that point times the valley
    %ratio
    max_value_im=subImaxPnts.*hess_bw.*sub_im*valleyRatio;
    dilate_kernel=true(minObjDim,minObjDim,minObjDim);
    subImax=imdilate(max_value_im,dilate_kernel);
    
   subImaxReg=subImaxPnts>subImax & sub_im>(2*thresh1);
    %make labelled mask
    hess_bwlabel=bwlabeln(hess_bw,6);
    %loop through labelled regions
    for iLabel=1:max(hess_bwlabel(:));
        %select specific region and the peaks in that region
        sub_hess_bw=hess_bwlabel==iLabel;
        subsubImax=subImaxReg & sub_hess_bw;
        
        %if more than one peak is found, watershed split them
        subsubcc=bwconncomp(subsubImax);
        if subsubcc.NumObjects>1
            maxBW=watershed(bwdist(subsubImax));
            hess_bw(maxBW==0 & subJm)=0;
        end
    end
end
%% after splitting, apply size filters
hess_bw=hess_bw.*sub_bw;
hess_bw=AreaFilter(hess_bw,minObjSize,[],6);
hess_bw=regionSplit(hess_bw,options);
% remove small objects in subImage
hess_bw=AreaFilter(hess_bw,minObjSize,[],6);


%% find centroids, display (off)
if show
centerIm(box(2):box(5),box(1):box(4),box(3):box(6))=Jc;
imagesc(sum(sub_im,3));
hold on
axis equal
temp=sum(Jc,3);
[y,x]=find(temp);
scatter(x,y,'black.');
hold off
pause(1);
end

%compile results of all sub images into final segmented mask.

im_bw_out(box(2):box(5),box(1):box(4),box(3):box(6))=...
    or(hess_bw,im_bw_out(box(2):box(5),box(1):box(4),box(3):box(6)));

end
%% 
