function [Vout, tformOut]=stackStabilization(V,searchRad,show,outputFlag)
% stackStabilization goes through a 3D image and aligns the points between
% consecutive frames. It starts at the middle frame and moves outward,
% finding peaks in each frame. It matches the peaks to the previous frame
% and calculates the affine transformation between the frames. It then
% accumulates the affine transforms from the middle so that a single
% transform can be applied to each slice


if nargin<4
    outputFlag=true;
end
if nargin<3
    show=0;
end
if nargin<2;
    searchRad=10;
    
end
if show
    outputFlag=1;
end
   
%
thresh=.01;
sliceLag=1;
V=normalizeRange(V);
imSize=size(V);
param.dim=2;
param.excessive=4;
param.quiet=1;
param.difficult=7.e2;
R = imref2d(size(V(:,:,1))) ;

Vtop=V(:,:,floor(imSize(3)/2+1):end);
Vbot=V(:,:,1:floor(imSize(3)/2));
Vout={Vbot,Vtop};
maxI=max(V(:));
%split image in half and do stabilization starting at the center.
for iHalf=1:2
    if iHalf==1
        Vtemp=flip(Vbot,3);
        
    else
        Vtemp=Vtop;
    end
    tformAll=[];
    VtempOut=Vtemp;
    VsmoothHistory=zeros(size(Vtemp));
    pointsHistory=cell(1,size(Vtemp,3));

    for iSlice=1:size(Vtemp,3)
                    pointsA=[];pointsB=[];

        imgBRaw=(Vtemp(:,:,iSlice))/maxI;
        imgBRaw(isnan(imgBRaw))=0;
        imgB=imgBRaw;
        if nnz(imgB)
        imgB=bpass(imgB,3,20);
        imgB=(imgB);
        imgB(imgB<thresh)=0;
        VsmoothHistory(:,:,iSlice)=imgB;
            
            [pointsB(:,2),pointsB(:,1)]=find(imregionalmax(imgB));
            pointsHistory{iSlice}=pointsB;
            
        if iSlice>sliceLag
            imgA= VsmoothHistory(:,:,iSlice-sliceLag);
            imgA=(imgA);
            imgA(imgA<thresh)=0;
            pointsA= pointsHistory{iSlice-sliceLag};

        else
            imgA=0;
        end
        else
            imgA=0;
            imgB=0;
        end
        
        if  nnz(imgB) && nnz(imgA) && iSlice>sliceLag

            trackInput=[pointsA  ones(size(pointsA(:,1))); ...
                pointsB 2*ones(size(pointsB(:,1)))];
            
            TrackOut=nan;
            for iSearch=0:searchRad
                if all(isnan(TrackOut(:)))
                    TrackOut=trackJN(trackInput,searchRad-iSearch,param);
                end
            end
            trackCounter=1;
            pointsAM=[];pointsBM=[];
            for i=1:TrackOut(end,end)
                trackItems=find(TrackOut(:,end)==i);
                if length(trackItems)==2
                    pointsAM(trackCounter,:)=TrackOut(trackItems(1),1:2);
                    pointsBM(trackCounter,:)=TrackOut(trackItems(2),1:2);
                    trackCounter=trackCounter+1;
                end
            end
            
            if length(pointsBM)>2
                
                %                 [tform, pointsBm, pointsAm] = estimateGeometricTransform(...
                %                     pointsBM, pointsAM, 'affine');
                
                tform=fitgeotrans( pointsBM, pointsAM, 'similarity');
%                 if det(tform.T)<.9 || det(tform.T)>1.1
%                     tform=affine2d(eye(3));
%                     
%                 end
                
            else
                tform=affine2d(eye(3));
            end
            
            %   tform.T(1:2,1:2)= tform.T(1:2,1:2)/det(tform.T);
            tform.T=tform.T*tformAll{iSlice-sliceLag}.T;
                        tformAll{iSlice}=tform;

            if outputFlag
            VtempOut(:,:,iSlice)=imwarp(Vtemp(:,:,iSlice),R,tform,'nearest',...
                'OutputView',R);
            if show
                imagesc(VtempOut(:,:,iSlice)-Vtemp(:,:,iSlice))
                hold on
                colormap jet
                if ~isempty(pointsAM)
                    scatter( pointsAM(:,1), pointsAM(:,2))
                    scatter( pointsBM(:,1), pointsBM(:,2))
                end
                hold off
                drawnow
            end
            else
                VtempOut(:,:,iSlice)=Vtemp(:,:,iSlice);
            end
        else
            tform=affine2d(eye(3));
            if iSlice>1
                tform.T=tform.T*tformAll{iSlice-1}.T;
            end
            tformAll{iSlice}=tform;
            
        end
    end
    
    if iHalf==1
        Vout=flip(VtempOut,3);
        tformOut=fliplr(tformAll);
    else
        Vout=cat(3,Vout,VtempOut);
        tformOut=[tformOut,tformAll];
    end
    
    
end


%
% [featuresA, pointsA] = extractFeatures(imgARaw, fliplr(pointsA));
% [featuresB, pointsB] = extractFeatures(imgBRaw, fliplr(pointsB));
%
% indexPairs = matchFeatures(featuresA, featuresB);
% pointsA = pointsA(indexPairs(:, 1), :);
% pointsB = cpointsB(indexPairs(:, 2), :);
% figure; showMatchedFeatures(imgA, imgB, pointsA, pointsB);
% legend('A', 'B');
%
%
