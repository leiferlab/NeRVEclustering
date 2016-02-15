function BWout=regionSplit(BW,options)
% recurrent function that splits up objects in a BW image if they are not
% sufficiently sphereical or have a long axis which is much longer than the
% other axes. 
minObjSize=60; % min object size
maxObjSize=350; % max object size
minSphericity=.83; % minimum sphericity for splitting.

% parse options to load fields
if nargin>=2
    Fnames=fieldnames(options);
    for i=1:length(Fnames)
        eval([Fnames{i} '= options.' Fnames{i} ';']);
    end
else
    options=struct;
end
options2=options;
options2.minSphericity=minSphericity*.9;
options2.maxObjSize=maxObjSize*1.5;

if isstruct(BW)
    cc=BW;
    BW=false(cc.ImageSize);
    BW(cell2mat(cc.PixelIdxList'))=true;
else
    cc=bwconncomp(BW,6);
end
stats=regionprops(cc,'BoundingBox');
% calculate blob parameters
perim_im=bwperim(BW,26);
vol=cellfun(@(x) length(x),cc.PixelIdxList);
perim=cellfun(@(x) sum(perim_im(x)),cc.PixelIdxList);
sphericity=((pi*36*vol.^2).^(1/3)./perim);
[D,V,Centroids]=regionPCA(cc);
D=(D).^(1/3);
longAxis=max(D,[],2);
Ddiff=D(:,3)-D(:,1);
aspectRatioCheck=Ddiff>max(6,D(:,1)*.7);

sizeCheck=vol>maxObjSize;
sphereCheck= sphericity<minSphericity;
%sphereCheck2=vol>maxObjSize*.8 && sphericity<minSphericity*
check=(aspectRatioCheck' | sphereCheck | sizeCheck) & ...
    vol>minObjSize & sphericity<1;



for iRegion=1:(cc.NumObjects)
    if check(iRegion)
        subIm=false(cc.ImageSize);
        subIm(cc.PixelIdxList{iRegion})=true;
        blank=subIm;
        watershedthresh=.4;
        waterCounter=0;
        newS=sphericity(iRegion);
        
        nObjs=1;
        [blankCrop,cropLookup]=imcrop3d(subIm,stats(iRegion).BoundingBox);
        Vsub=squeeze(V(iRegion,:,:));
        Jw=ones(size(blankCrop));
        while((all(Jw(:))) || ~sum(~Jw(blankCrop))) && watershedthresh>.2
            
            
            
            
            if ~waterCounter
                if aspectRatioCheck(iRegion)
                    % still try to seperate objects that are flat, but
                    % doing this to make BWdist not acccount for the thin
                    % dimension as much;
                    [~,dmin]=max(abs(Vsub),[],2);
                    correctionScale=1./D(iRegion,:);
                    %ignore flat objects
                    if all(isfinite(correctionScale))
                    correctionScale=correctionScale/min(correctionScale);
                    correctionScale=correctionScale(dmin);
                    correctionScale(dmin==1)=5*correctionScale(dmin==1);
                    correctionScale(dmin==2)=correctionScale(dmin==2);
                    correctionScale=round(correctionScale);
                    Jd=-bwdist_jn(~blankCrop,correctionScale);
                    else
                        Jd=zeros(size(blankCrop));
                    end
                else
                    Jd=-bwdist_jn(~blankCrop);
                end
                waterCounter=true;
            end
            
            
            %Jd=smooth3(Jd,'gaussian',5,2);
            Jdh=imhmin(Jd,watershedthresh);
            Jdh(~blankCrop)=Inf;
            subsubcc=bwconncomp(imregionalmin(Jdh));
            
            if subsubcc.NumObjects>1
                
                Jw=watershed(Jdh);
                blank2=blankCrop;
                blank2(blankCrop)=~~Jw(blankCrop);
                blankPerim=bwperim(blank2,26);
                Jw_cc=bwconncomp(blank2,6);
                perims=cellfun(@(x) sum(blankPerim(x)),Jw_cc.PixelIdxList);
                vols=cellfun(@(x) length(x),Jw_cc.PixelIdxList);
                newS=min((pi*36*vols.^2).^(1/3)./perims);
                nObjs=subsubcc.NumObjects;
            else
                if  subsubcc.NumObjects==1
                    blank2=blankCrop;
                end
            end
            
            
            watershedthresh=watershedthresh-.1;
            
        end
        %seperate regions and then run regionSplit again on seperated
        %regions
        if nObjs>1
            subIm(blank)=blank2(blankCrop(:));
            subIm=regionSplit(subIm,options);
                    BW(cc.PixelIdxList{iRegion})=subIm(blank);
        elseif nObjs==1 && ( aspectRatioCheck(iRegion) || newS<minSphericity ||vol(iRegion)>maxObjSize*1.5)
            [blanky,blankx,blankz]=ind2sub(cc.ImageSize,cc.PixelIdxList{iRegion});
            blankCentered=bsxfun(@minus,[blankx blanky blankz],Centroids(iRegion,:));
            %  blankCentered=blankCentered.^2
            blankProj=abs(blankCentered*Vsub([2 1 3],end));
            subIm(cc.PixelIdxList{iRegion}(blankProj<.6))=false;
                        subIm=regionSplit(subIm,options2);
        BW(cc.PixelIdxList{iRegion})=subIm(blank);
        end
        
    end
    
end
BW=AreaFilter(BW,minObjSize,[],6);
BWout=BW;

