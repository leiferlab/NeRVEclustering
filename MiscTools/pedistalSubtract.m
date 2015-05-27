function [ outputImage,pedistal ] = pedistalSubtract( inputImage,border)
%pedistalSubtract takes an image and subtracts the median of the values
%around the border as a form of pedestal subtraction for images. The border
%default size is 3.

% passing a "border" argument of "zn" (z1,z2,z3 ...) fits a polynomials to
% the order listed and subtracts that from the stack


% check the first input
if nargin == 1
    border =3;
    borderString = 'z0';
end

% check the second input
if nargin == 2;
    if ischar(border)
       borderString = border;
       border = 3;
       
    else
        borderString = 'z0';
    end
end


% find the background of the image as the border voxels
pedMaskf=false(size(inputImage));
% pull out the voxels that are around the borders
pedMaskf(1:border,:,:)=true;
pedMaskf(:,1:border,:)=true;
pedMaskf(end-border+1:end,:,:)=true;
pedMaskf(:,end-border+1:end,:)=true;
% suppress any voxels that are equal to the minimum, for instance if there
% are big swatches of zeros from some rotation
pedMaskf(inputImage==min(inputImage(:)))=false;
% find the pedestal for each plane
maskedInputImage = double(inputImage);
maskedInputImage(~(pedMaskf)) = nan;
for iiPlane = 1:size(inputImage,3);
    tempPlane = maskedInputImage(:,:,iiPlane);
   pedestal(iiPlane) = nanmedian(tempPlane(:)); 
end

% find what order polynomial to fit
polyOrder = str2double(borderString(2:end));
if isnan(polyOrder);
    warning('pedestalSubtract:wrongType','''Border'' input is not of the correct form');
    polyOrder = 0;
end

pedSelect=find(~isnan(pedestal));
% fit to a polynomial
pedestalCoefs = polyfit(pedSelect,pedestal(pedSelect),polyOrder);
% evalutate the polynomial
pedestalVals = polyval(pedestalCoefs,1:size(inputImage,3));
% shape it along the third dimension
pedestalVals = reshape(pedestalVals,1,1,[]);
% pedistal=median(inputImage(pedMaskf));
% subtract this value from each plane
outputImage = bsxfun(@minus,inputImage,pedestalVals);
outputImage(outputImage<0) = 0;

% save the pedestal values
if nargout>1
   pedistal = pedestalVals; 
end

end

