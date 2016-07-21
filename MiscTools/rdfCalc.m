function [counts2,radialDistribNew,rEvalList,thetaList]=rdfCalc(inputImage,params);
% Input parameters are
% params is a structure with the following fields:
% wrapEdge is 'none','TB','LR', or 'both' which defines whether to consider
%   the image to have periodic boundary conditions on any of the sides
% radialStepSize is the binning size for RDF evaluation in pixels as in,
%   the radius for the distance between two pixels is rounded to the nearest
%   radialStepSize
% angularStepSize is the binning size for the angular component of the RDF
% longDistance, should be set to FALSE
% calcBruteForce, should be set to FALSE unless binary image
% calcConvMethod, should be set to TRUE
% displayPlots, should be set to FALSE

%% simulate image (must be even dimensions?)
if nargin<1
    imageSize = [15,15];
    inputImage = zeros(imageSize);
    nPts = 10;
    randomPts = [];
    randomPts(:,1) = round(1+(imageSize(1)-1)*rand(nPts,1));
    randomPts(:,2) = round(1+(imageSize(2)-1)*rand(nPts,1));
    %
    % random positions
    indices = sub2ind(size(inputImage),randomPts(:,1),randomPts(:,2));
    inputImage(indices) = 1;
    
    %     % gridded positions
    %     inputImage = zeros(size(inputImage));
    %     inputImage(2:25:end,2:25:end) = 1;
    
    % % edges only
    % % inputImage = zeros(size(inputImage));
    % inputImage(:,1) = 1;
    % inputImage(:,end) = 1;
    % inputImage(1,:) = 1;
    % inputImage(end,:) = 1;
    
    
    % % corners only
    % inputImage = zeros(size(inputImage));
    % inputImage(1,1) = 1;
    % inputImage(1,end) = 1;
    % inputImage(end,1) = 1;
    % inputImage(end) = 1;
    %
    
    % small disks
    %     inputImage = imdilate(inputImage,strel('disk',2,0));
    
end

% coordinates of bright pixels
randomPts2 = nan(nnz(inputImage(:)),2);
[randomPts2(:,1),randomPts2(:,2)] = find(inputImage);
randomPts = randomPts2;
nPts = numel(randomPts)/2;

%%
if nargin<2
    % default parameters
    params.wrapEdge = 'none';
    params.radialStepSize = 0.05;
    params.angularStepSize = pi/2;
    params.longDistance = false;
    params.calcBruteForce = true;
    params.calcConvMethod = true;
    params.displayPlots = true;
end

%% process image
if params.displayPlots
    clf;
end
imageSize = size(inputImage);

% generate grid
xSize = size(inputImage,2);
xList = -1*xSize:xSize;
ySize = size(inputImage,1);
yList = -1*ySize:ySize;
[xList,yList] = meshgrid(xList,yList);

if params.longDistance;
    % biggest radius to evaluate is the short edge length
    biggestREval = sqrt(2)*min(imageSize)/2;
else
    switch params.wrapEdge
        case 'both'
            % biggest radius is the shortest side (matches brute force)
            biggestREval = min(imageSize)/2;
        case 'LR'
            biggestREval = min(imageSize(1),imageSize(2)/2);
        case 'TB'
            biggestREval = min(imageSize(1)/2,imageSize(2));
        otherwise % including false
            biggestREval = min(imageSize)-1;
    end
    
end

% list of radii to bin into
rEvalList = 0.5:params.radialStepSize:biggestREval;

% list of angles to bin into
thetaList = -pi/2:params.angularStepSize:pi/2;
if numel(thetaList)<2
    thetaList = pi/4;
    params.angularStepSize = pi/2;
end

if params.calcConvMethod
    % calculate radii for pixel shifts
    rTemp = sqrt((xList).^2+(yList).^2);
    rMask = nan(size(inputImage)*2+1);
    
    % make a mask of distance labels
    for iREval = 1:numel(rEvalList)-1
        rMask(and(rTemp>rEvalList(iREval),rTemp<=rEvalList(iREval+1))) = iREval;
    end
    
    % calculate angles for orientational RDF
    
    % absolute value averages chiralities
%     thetaTemp = abs(atan(yList./xList));
    thetaTemp = atan(yList./xList);
    
    thetaMask = nan(size(inputImage)*2+1);
    
    % make a mask of angle labels
    for iTheta = 1:numel(thetaList)
        if iTheta == numel(thetaList) % include pi/2 in last bin
            thetaMask(and(thetaTemp>=(thetaList(iTheta)-params.angularStepSize/2),...
            thetaTemp<=thetaList(iTheta)+params.angularStepSize/2)) = iTheta;
        else  
        thetaMask(and(thetaTemp>=(thetaList(iTheta)-params.angularStepSize/2),...
            thetaTemp<thetaList(iTheta)+params.angularStepSize/2)) = iTheta;
        end
    end
    % replace nans with a 0
    thetaMask(isnan(thetaMask))=0;
    
    % depending on the choice, replicate the input matrix to deal with edge
    % effects
    % generate a mask for proper normalization
    inputMask = zeros(3*size(inputImage));
    inputMask(ySize+1:2*ySize,xSize+1:2*xSize)=1*sum(inputImage(:));
    
    switch params.wrapEdge
        case 'LR' % left and right
            % replicate image to be the same size
            inputImgRep = zeros(3*size(inputImage,1),3*size(inputImage,2),2);
            % singleton
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                size(inputImage,2):(2*size(inputImage,2))-1,1) = inputImage;
            % rep horizontal
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                round(size(inputImage,2)/2):round((5*size(inputImage,2)/2)-1),2) = ...
                repmat(inputImage,[1,2]);
            % image to display
            inputImgRepShow = inputImgRep(:,:,2);
        case 'TB' % top and bottom
            % replicate image to be the same size
            inputImgRep = zeros(3*size(inputImage,1),3*size(inputImage,2),2);
            % singleton
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                size(inputImage,2):(2*size(inputImage,2))-1,1) = inputImage;
            % rep vertical
            inputImgRep(round(size(inputImage,1)/2):round((5*size(inputImage,1)/2-1)),...
                size(inputImage,2):(2*size(inputImage,2))-1,2) = ...
                repmat(inputImage,[2,1]);
            % image to display
            inputImgRepShow = inputImgRep(:,:,2);
        case 'both' % periodic in 2D
            % replicate image to be the same size
            inputImgRep = zeros(3*size(inputImage,1),3*size(inputImage,2),4);
            % singleton
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                size(inputImage,2):(2*size(inputImage,2))-1,1) = inputImage;
            % rep vertical
            inputImgRep(round(size(inputImage,1)/2):round((5*size(inputImage,1)/2-1)),...
                size(inputImage,2):(2*size(inputImage,2))-1,2) = ...
                repmat(inputImage,[2,1]);
            % rep horizontal
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                round(size(inputImage,2)/2):round((5*size(inputImage,2)/2))-1,3) = ...
                repmat(inputImage,[1,2]);
            % rep both vert and horizontal
            inputImgRep(round(size(inputImage,1)/2+1):round((5*size(inputImage,1)/2)),...
                round(size(inputImage,2)/2+1):round((5*size(inputImage,2)/2)),4) = ...
                repmat(inputImage,[2,2]);
            
            % image to display
            inputImgRepShow = inputImgRep(:,:,4);
        otherwise % including the false case
            % replicate image to be the same size
            inputImgRep = zeros(3*size(inputImage,1),3*size(inputImage,2),1);
            % singleton
            inputImgRep(size(inputImage,1):(2*size(inputImage,1)-1),...
                size(inputImage,2):(2*size(inputImage,2))-1,1) = inputImage;
            % image to display
            inputImgRepShow = inputImgRep(:,:,1);
    end
    
    if params.displayPlots;
        % display input image, zero padded and replicated
        subplot(1,2,1);
        imshow(inputImgRepShow,[])
        title('original image (replicated according to periodicity');
    end;
    
    % preallocate arrays
    convImg = zeros(size(inputImgRep));
    radialDistrib = nan(numel(rEvalList),5,numel(thetaList));
    
    % area density for future normalization
    areaDensity = sum(inputImage(:))./numel(inputImage);
    
    %
    yCenter = size(inputImage,1)+1;
    xCenter = size(inputImage,2)+1;
    % loop over each radius in rEvalList and calculate p(r) using convolution
    for iREval = 1:numel(rEvalList)-1
        for iTheta = 1:numel(thetaList)
            currentRadius = ceil(rEvalList(iREval));
            %     determine which pixels will be used for this radius
            ringImg = double(...
                rMask(yCenter-currentRadius:yCenter+currentRadius,...
                xCenter-currentRadius:xCenter+currentRadius)==iREval);
            % determine which pixels will be used for this angle
            ringImg = ringImg.*double(...
                thetaMask(yCenter-currentRadius:yCenter+currentRadius,...
                xCenter-currentRadius:xCenter+currentRadius)==iTheta);
            
            %         imshow(ringImg,[],'InitialMagnification','fit');
            %         figure(gcf);
            %         pause(0.5)
            %         ringImg = double(rMask==iREval);
            %         imshow(ringImg,[],'InitialMagnification','fit');
            %         figure(gcf);
            %         pause(0.5)
            ringTot = nnz(ringImg(:));
            
            % pick up values in circles using convolution
            for kkRep = 1:size(inputImgRep,3);
                convImg = conv2(inputImgRep(:,:,kkRep),ringImg,'same');
                
                % weight convolution by initial pixel values
                prodImg = convImg.*inputImgRep(:,:,kkRep);
                
                
                % radial distribution is a ratio of sums, this only
                % calculates one part of that ratio
                radialDistrib(iREval,kkRep,iTheta) = sum(prodImg(:));
            end
        end
    end
    
    % deal with overlaps for periodicity
    switch params.wrapEdge
        case {'TB','LR'}
            radialDistribNew = radialDistrib(:,2,:)-radialDistrib(:,1,:);
        case 'both'
            radialDistribNew = radialDistrib(:,1,:)+radialDistrib(:,4,:)-radialDistrib(:,2,:)-radialDistrib(:,3,:);
        otherwise
            radialDistribNew = radialDistrib(:,1,:);
    end
    radialDistribNew = squeeze(radialDistribNew);
    % the initial point of zero overlap is the pointyness of the original
    % image
    radialDistribNew(1,:) = sum(inputImage(:).^2);
    
    if params.longDistance
        warning('rdfCalculator:longDistance', 'Longer distance calculations should be avoided, they do not seem to match the brute force method exactly');
%         % enumerate overlaps for correcting image convolution
%         counter = 0;
%         resultsA = [];
%         for ii = 0:imageSize(1) % distance 1
%             for jj = 0:imageSize(1) % distance 2
%                 for kk = 0:imageSize(2) % distance 3
%                     for ll = 0:imageSize(2) % distance 4
%                         ij = ii+jj;
%                         kl = kk+ll;
%                         % offsets add to proper size
%                         if not(or(and(ij==imageSize(1),kk==ll),and(ii==jj,kl==imageSize(2))));
%                             continue
%                         end
%                         
%                         % one step not too long
%                         if sqrt(ii^2+kk^2)>(min(imageSize)/2)*sqrt(2);
%                             continue
%                         end
%                         
%                         % other step not too long
%                         if sqrt(jj^2+ll^2)>(min(imageSize)/2)*sqrt(2);
%                             continue
%                         end
%                         
%                         counter = counter+1;
%                         resultsA(counter,:) = [ii,jj,kk,ll,sqrt(ii^2+kk^2),sqrt(jj^2+ll^2)];
%                         
%                     end
%                 end
%             end
%         end
%         
%         % calculate the values to subtract off
%         distancesToUse = resultsA(:,5:6);
%         distancesToCalc = max(distancesToUse,[],2);
%         distancesToCalc = unique(distancesToCalc);
%         distancesToCalc = distancesToCalc(1:end-1);
%         nDists = numel(distancesToCalc);
%         
%         radialDistrib2 = [];
%         for iSubtract = 1:nDists;
%             % find index in main list
%             try
%                 mainInd(iSubtract) = find(rEvalList>distancesToCalc(iSubtract),1,'first')-1;
%                 
%                 % find companion length
%                 resultsB = cat(1,...
%                     resultsA(resultsA(:,5)==distancesToCalc(iSubtract),5:6),...
%                     resultsA(resultsA(:,6)==distancesToCalc(iSubtract),5:6));
%                 resultsB = min(resultsB,[],2);
%                 resultsB = unique(resultsB);
%                 
%                 if numel(resultsB)~=1
%                     % very special case, treat separately
%                     variableType = 1;
%                 elseif (resultsB-floor(resultsB))<1e-12;
%                     % integers, remove as special case
%                     variableType = 2;
%                 elseif resultsB == distancesToCalc(iSubtract);
%                     % doubled values, remove half of calculation
%                     variableType = 3;
%                 else
%                     % everything else, remove entire calculation
%                     variableType = 4;
%                 end
%                 
%                 % create ring image for convolution
%                 ringImg = double(rTemp==distancesToCalc(iSubtract));
%                 
%                 %     pause(0.1);
%                 %             ringTot = nnz(ringImg(:));
%                 % pick up values in circles using convolution
%                 for kkRep = 1:size(inputImgRep,3);
%                     convImg(:,:,kkRep) = conv2(inputImgRep(:,:,kkRep),ringImg,'same');
%                     
%                     %     imshow(inputImgRep(:,:,5),[])
%                     % weight convolution by initial pixel values
%                     prodImg = convImg(:,:,kkRep).*inputImgRep(:,:,kkRep);
%                     
%                     % NEEDS WORK
%                     % normalization needs to take into consideration the edge effects and
%                     % limited area that contributes to each pixel
%                     %                 normImg = conv2(inputMask,ringImg,'same').*inputImgRep(:,:,kkRep);
%                     
%                     % radial distribution is a ratio of sums
%                     radialDistrib2(iSubtract,kkRep) = sum(prodImg(:));%./(sum(normImg(:)))./areaDensity;
%                     
%                 end
%                 
%                 % deal with overlaps for periodicity
%                 switch params.wrapEdge
%                     case {'TB','LR'}
%                         radialDistribNew2 = radialDistrib2(:,2)-radialDistrib2(:,1);
%                     case 'both'
%                         radialDistribNew2 = radialDistrib2(:,1)+radialDistrib2(:,4)-radialDistrib2(:,2)-radialDistrib2(:,3);
%                     otherwise
%                         radialDistribNew2 = radialDistrib2(:,1);
%                 end
%                 
%                 
%                 % phenomenomological
%                 radialDistribNew2 = radialDistribNew2/3;
%                 
%                 % phenomenomological correction remove from correct bin
%                 switch variableType
%                     case {4}
%                         radialDistribNew(mainInd(iSubtract)) = ...
%                             radialDistrib(mainInd(iSubtract)) - ...
%                             radialDistribNew2(iSubtract)/8;
%                     case {1}
%                         radialDistribNew(mainInd(iSubtract)) = ...
%                             radialDistrib(mainInd(iSubtract)) - ...
%                             radialDistribNew2(iSubtract)/32;
%                     case {2}
%                         radialDistribNew(mainInd(iSubtract)) = radialDistrib(mainInd(iSubtract)) - ...
%                             radialDistribNew2(iSubtract)/16;
%                     case {3}
%                         radialDistribNew(mainInd(iSubtract)) = radialDistribNew(mainInd(iSubtract)) - ...
%                             radialDistribNew2(iSubtract);
%                     otherwise
%                         radialDistribNew(mainInd(iSubtract)) = 0;
%                 end
%             catch ME
%                 distancesToCalc = distancesToCalc(1:iSubtract-1);
%                 break
%             end
%         end
%     else
%         
    end
else
    radialDistribNew = [];
end % calculate convolution method

%%  brute force calculate rdf (dealing with periodic boundary condtions)
if params.calcBruteForce
    distances2 = pdist2(randomPts,randomPts);
    periodicOffset1 = zeros(nPts,2);
    periodicOffset1(:,1) = 1;
    periodicOffset2 = -1*(periodicOffset1-1);
    
    % calculate distances using proper periodicity
    switch params.wrapEdge
        case 'LR'
            distances4 = pdist2(randomPts,randomPts-imageSize(2)*periodicOffset2);
            distances7 = pdist2(randomPts,randomPts+imageSize(2)*periodicOffset2);
            distances = cat(3,distances2,distances4,distances7);
        case 'TB'
            distances3 = pdist2(randomPts,randomPts-imageSize(1)*periodicOffset1);
            distances6 = pdist2(randomPts,randomPts+imageSize(1)*periodicOffset1);
            distances = cat(3,distances2,distances3,distances6);
        case 'both'
            distances3 = pdist2(randomPts,randomPts-imageSize(1)*periodicOffset1);
            distances4 = pdist2(randomPts,randomPts-imageSize(2)*periodicOffset2);
            distances5 = pdist2(randomPts,randomPts-imageSize(1)*periodicOffset1-imageSize(2)*periodicOffset2);
            distances6 = pdist2(randomPts,randomPts+imageSize(1)*periodicOffset1);
            distances7 = pdist2(randomPts,randomPts+imageSize(2)*periodicOffset2);
            distances8 = pdist2(randomPts,randomPts+imageSize(1)*periodicOffset1+imageSize(2)*periodicOffset2);
            distances9 = pdist2(randomPts,randomPts-imageSize(1)*periodicOffset1+imageSize(2)*periodicOffset2);
            distances10 = pdist2(randomPts,randomPts+imageSize(1)*periodicOffset1-imageSize(2)*periodicOffset2);
            distances = cat(3,distances2,distances3,distances4,distances5,distances6,...
                distances7,distances8,distances9,distances10);
        otherwise
            distances = distances2;
    end
    
    [minDist,minInd] = min(distances,[],3);
    
    counts2 = histc(minDist(:),rEvalList);
    radialDistrib2Area = pi*((rEvalList(2:end).^2-...
        (rEvalList(1:end-1).^2)));
    radialDistrib2Area = [rEvalList(1).^2,radialDistrib2Area];
    totalPts = nnz(inputImage);
    
else
    counts2 = [];
end
% normalized result
if params.displayPlots;
    subplot(1,2,2);
    cla;
    hold on;
    if params.calcBruteForce
        scatter(rEvalList,counts2','ro');
    end
    if params.calcConvMethod
        scatter(rEvalList,radialDistribNew,'kx');
    end
    grid on;
    title('radial distribution function (counts)');
end

end


