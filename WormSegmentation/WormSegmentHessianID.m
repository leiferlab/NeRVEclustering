function wormBW2=WormSegmentHessianID(worm,options)
%WORMSEGMENTHESSIANID takes a 3D matrix representation of an image stack
%and a structure of parameter fields and outputs a matrix representation of
%the binary image stack with significant objects separated.

%NOTE: requires many functions from 3dbrain directory to function

%% 

%Initialize default parameters; all of these are also fields in options
thresh1 = .03;              %initial threshold
hthresh = -.0001;           %threshold for trace of hessian
minObjSize = 500;           %min object size
maxObjSize = 5000;          %max object size
watershedFilter = 1;        %watershed filter for object shapes
filterSize = [20,20,10];    %bandpass filter size low f
noise = 1;                  %bandpass filter hi f
pad = 4;                    %pad to take around each sub blob
maxSplit = 1;               %split objects using regional maxima
minSphericity = .55;        %minimum sphericity for splitting

%parse options to load fields
if nargin == 2                          %if "options" specified...
    Fnames = fieldnames(options);
    
    for i = 1:length(Fnames)            %for each parameter in "options"...
        eval([Fnames{i} '= options.' Fnames{i} ';']);
    end
    
end

%% 

%Create array of dimensions of the image stack
imsize = size(worm);

imsize = imsize([2,1,3]); %                                                       <-- what is the point of this statement??

%Construct pedistal (perimeter of image space)
pedMask = false(imsize);
pedMask(1:3,:,:)       = true;
pedMask(:,1:3,:)       = true;
pedMask(end-2:end,:,:) = true;
pedMask(:,end-2:end,:) = true;

%Subtract pedistal as a background value
pedistal = median(worm(pedMask));
worm = worm - pedistal;
worm(worm < 0) = 0;

%Pass through bandpass filter, normalize and threshold intensities
wormtop = bpass3_jn(worm,noise,filterSize);
wormtop = normalizeRange(wormtop);
wormBW = wormtop > thresh1;

%% 

%Find connected components, ignore those smaller than specified minimum
cc = bwconncomp(wormBW,6);
blobSizes = cellfun(@(x) length(x), cc.PixelIdxList);
cc.PixelIdxList(blobSizes < minObjSize) = [];
cc.NumObjects=sum(blobSizes >= minObjSize);
blobStats = regionprops(cc,'Area','BoundingBox','Centroid');

%Use hessian to find nuclei in objects
wormBW2 = zeros(size(wormBW));

for iblob = 1:cc.NumObjects;            %for each connected component...
    
    %Crop out current object with pad
    BB = floor(blobStats(iblob).BoundingBox);
    BB(1:3) = BB(1:3) - [pad,pad,pad];
    BB(4:end) = BB(4:end) + BB(1:3) + [2*pad,2*pad,2*pad];
    
    %Don't crop area outside original image
    BB(BB < 1) = 1;
    BB([false,false,false,bsxfun(@ge,BB(4:6),imsize)]) = imsize((bsxfun(@ge,BB(4:6),imsize)));
    
    %Create and normalize sub-image containing current object
    subBW = wormBW(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    subIm = wormtop(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6));
    subIm = normalizeRange(subIm);
    
    %Smooth image and calculate hessian and eigenvalues for segmentation
    subIm = smooth3(subIm,'gaussian',5,3);
    H = hessianMatrix(subIm,3);
    Heig = hessianEig(H);
    Htrace = max(Heig,[],4);
    Htrace(isnan(Htrace)) = 0;
    Jm = Htrace < hthresh;
    Jm = xyzConvHull(Jm,3);
    
    %Watershed filter shapes
    if watershedFilter
        Jd = -bwdist(~Jm);
        Jd = imhmin(Jd,watershedFilter);
        Jd(~Jm) = Inf;
        Jw = watershed(Jd);
        Jm = Jm.*(Jw > 0);
    end
    
    Jm = imerode(Jm,true(2,2,2));
    Jm = imclearborder(Jm);
    
    
    if maxSplit    %Watershed splitting based on local maxima locations

        %Cuberoot minimum size parameter for minimum object dimension
        minObjDim = round(minObjSize^(1/3));
        
        %Create binary image of regional maxima
        subImax = imregionalmax(subIm);
        subImax = imdilate(subImax.*Jm,true(minObjDim,minObjDim,minObjDim));
        
        %Label each connected component and adjust binary image accordingly
        JmLabel = bwlabeln(Jm,6);
        for iLabel = 1:max(JmLabel(:));     %for each labeled component...
            subJm = JmLabel==iLabel;
            subsubImax = subImax & subJm;
            maxBW = watershed(bwdist(subsubImax));
            Jm(maxBW==0 & subJm) = 0;
        end
        
    end
    
    %Find connected components in newly-watershed binary image
    subCC = bwconncomp(Jm>0,6);
    for isubBlob = 1:subCC.NumObjects %for each sub-connected-component...
        
        %If current object meets specified minimum size criteria...
        if length(subCC.PixelIdxList{isubBlob}) > minObjSize
           
            %Check volume and sphericity before further watershed
            blank = false(size(subBW));
            blank(subCC.PixelIdxList{isubBlob}) = true;
            subblobP = bwperim(blank,26);                     %perimeter
            subblobP = sum(subblobP(:));    
            
            subblobV = sum(blank(:));                         %volume
            subblobS = (pi*36*subblobV.^2)^(1/3) / subblobP;  %sphericity
            
            %If current object fails either size or sphericity criteria...
            if subblobV > maxObjSize || subblobS < minSphericity
                Jw=ones(size(Jd));
                watershedthresh=.7;
                
                while ((all(Jw(:))) || ~sum(~Jw(blank))) && watershedthresh > .4
                    Jd = -bwdist(~blank);
                    Jd = imhmin(Jd,watershedthresh);
                    Jd(~Jm) = Inf;
                    Jw = watershed(Jd);
                    watershedthresh = watershedthresh - .1;
                end
                
                Jm(blank) = ~~Jw(blank);
            end
        end
        
    end
    
    % remove small objects in subImage
    Jm = AreaFilter(Jm,minObjSize, maxObjSize,6);

    %compile results of all sub images into final segmented mask.
    Jm = imclearborder(Jm);
    wormBW2(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6)) = or(Jm,wormBW2(BB(2):BB(5),BB(1):BB(4),BB(3):BB(6)));
    
end


end