function [Vout, tformOut]=stackStabilization(V,searchRad,show)
% stackStabilization goes through a 3D image and aligns the points to the
if nargin<3
    show=0;
end
if nargin<2;
    searchRad=10;
    
end
thresh=.01;
sliceLag=1;
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
for iHalf=1:2
    if iHalf==1
        Vtemp=flip(Vbot,3);
        
    else
        Vtemp=Vtop;
    end
    tformAll=[];
    VtempOut=Vtemp;
    
    for iSlice=1:size(Vtemp,3)
        
        imgBRaw=(Vtemp(:,:,iSlice))/maxI;
        imgBRaw(isnan(imgBRaw))=0;
          imgB=imgBRaw;
            imgB=bpass(imgB,3,20);
            imgB=(imgB);
            imgB(imgB<thresh)=0;
            if iSlice>sliceLag
            imgARaw=(Vtemp(:,:,iSlice-sliceLag))/maxI;
            imgA=imgARaw;
            imgA=bpass(imgA,3,20);
            imgA=(imgA);
            imgA(imgA<thresh)=0;
            else
                imgA=0;
            end      
        
        if  nnz(imgB) && nnz(imgA) && iSlice>sliceLag
            %%
     

            
            % pointsA = detectFASTFeatures(imgA, 'MinContrast', ptThresh);
            % pointsB = detectFASTFeatures(imgB, 'MinContrast', ptThresh);
            pointsA=[];pointsB=[];
            [pointsA(:,2),pointsA(:,1)]=find(imregionalmax(imgA));
            [pointsB(:,2),pointsB(:,1)]=find(imregionalmax(imgB));
            
            % if show
            %     %%
            % subplot(1,2,1); imagesc(imgA); hold on;
            % scatter(pointsA(:,1),pointsA(:,2));
            % title('Corners in A');
            % hold off
            %
            %
            % subplot(1,2,2); imagesc(imgB); hold on;
            % scatter(pointsB(:,1),pointsB(:,2));
            % title('Corners in B');
            % hold off
            % pause(.2)
            % end
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

tform=fitgeotrans( pointsBM, pointsAM, 'affine');
            else
                tform=affine2d(eye(3));
            end
          %  [det(tform.T) length(pointsBM)]
            tform.T=tform.T*tformAll{iSlice-sliceLag}.T;
            VtempOut(:,:,iSlice)=imwarp(Vtemp(:,:,iSlice),R,tform,...
                'OutputView',R);
            tformAll{iSlice}=tform;
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
