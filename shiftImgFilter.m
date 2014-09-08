function [img,tformOut] = shiftImgFilter(img,filterID,tformIn)
% INPUTS
% img is an m x n x 1 array containing pixel values
% filterID is the integer ID of the filter used to take the image
% tformIn is an optional parameter to pass in the the tform to be used from
% image registration.
% 
% OUTPUTS
% img is an m x n x 1 array shifted to a common location for all filters
% tformOut is the tform transform used to shift the image
%
% USAGE
% [img,tformOut] = shiftImgFilter(img,filterID);
% [img] = shiftImgFilter(img,[],tformIn);

% settings for various filters
if nargin<3 % no default tform to use
    switch filterID
        case 1
            x = -3.25;
            y = 1.30;
            scale = 1;
        case 2
            x = -1.54;
            y = 1.23;
            scale = 1;
        case 3 % considered the 'standard' filter
            x = 0;
            y = 0;
            scale = 1;
        case 4 
            x = -1.45;
            y = -0.66;
            scale = 1;
        case 5 
            % this is the only one that requires scaling and it is not
            % quite clear to BPB how to check and make sure it is working
            % well.
            x = 1.79;
            y = -0.44;
            scale = 1.0031667;
            % should the scale be 1/scale? x/y are -deltaX, -deltaY
        case 6
            x = -0.77;
            y = 0.15;
            scale = 1;
        otherwise
            x = 0;
            y = 0;
            scale = 1;
    end
    
    tformOut = maketform('affine',[[scale,0,0];[0,scale,0];[x,y,1]]);
else
    tformOut = tformIn;
    
end

% shift/transform the image but also change the coordinate system
img = imtransform(img,tformOut,...
    'XData',[1, size(img,2)],...
   'YData',[1, size(img,1)]);
