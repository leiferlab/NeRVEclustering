function res = bpass3_jn(image_array,lnoise,lobject,threshold, isForcedDouble)
% Implements a real-space bandpass filter that suppresses pixel noise and long-wavelength image variations while retaining information of a characteristic size.
%
%
% NAME:
%               bpass
% PURPOSE:
%               Implements a real-space bandpass filter that suppresses 
%               pixel noise and long-wavelength image variations while 
%               retaining information of a characteristic size.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = bpass( image_array, lnoise, lobject )
% INPUTS:
%               image:  The two-dimensional array to be filtered.
%               lnoise: Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%               lobject: (optional) Integer length in pixels somewhat 
%                       larger than a typical object. Can also be set to 
%                       0 or false, in which case only the lowpass 
%                       "blurring" operation defined by lnoise is done,
%                       without the background subtraction defined by
%                       lobject.  Defaults to false.
%               threshold: (optional) By default, after the convolution,
%                       any negative pixels are reset to 0.  Threshold
%                       changes the threshhold for setting pixels to
%                       0.  Positive values may be useful for removing
%                       stray noise or small particles.  Alternatively, can
%                       be set to -Inf so that no threshholding is
%                       performed at all.
%
% OUTPUTS:
%               res:    filtered image.f
% PROCEDURE:
%               simple convolution yields spatial bandpass filtering.
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.

% MODIFICATION HISTORY:
%               Written by David G. Grier, The University of Chicago, 2/93.
%
%               Greatly revised version DGG 5/95.
%
%               Added /field keyword JCC 12/95.
% 
%               Memory optimizations and fixed normalization, DGG 8/99.
%               Converted to Matlab by D.Blair 4/2004-ish
%
%               Fixed some bugs with conv2 to make sure the edges are
%               removed D.B. 6/05
%
%               Removed inadvertent image shift ERD 6/05
% 
%               Added threshold to output.  Now sets all pixels with
%               negative values equal to zero.  Gets rid of ringing which
%               was destroying sub-pixel accuracy, unless window size in
%               cntrd was picked perfectly.  Now centrd gets sub-pixel
%               accuracy much more robustly ERD 8/24/05
%
%               Refactored for clarity and converted all convolutions to
%               use column vector kernels for speed.  Running on my 
%               macbook, the old version took ~1.3 seconds to do
%               bpass(image_array,1,19) on a 1024 x 1024 image; this
%               version takes roughly half that. JWM 6/07
%
%               20100710 BPB. Refactored for speed and updated matlab algorthms, 
%               no longer faster to separate dimensions of kernels.
%
%               20100712 BPB. Added single/double switch. It is not yet an
%               optional keyword that can be anywhere, it may only come as
%               the fifth option.
%
%       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
%       David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  


if nargin < 3, lobject = false; end
if nargin < 4, threshold = 0; end
if nargin < 5, isForcedDouble = true; end
if numel(lobject)==1, lobject=repmat(lobject,1,3);end
normalize = @(x) x/sum(x(:));

if isForcedDouble
    image_array = double(image_array);
else
    imageType = class(image_array);
    if not(strcmp(imageType,{'single','double'}))
    image_array = single(image_array);
    end
end


 
% JWM: Do a 2D convolution with the kernels in two steps each.  It is
% possible to do the convolution in only one step per kernel with 
%
  % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
  % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
  % filtered = conv2(image_array,...
  %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
  %   'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.

% gconv = conv2(image_array',gaussian_kernel','same');
% gconv = conv2(gconv',gaussian_kernel','same');
% 
% if lobject
%   bconv = conv2(image_array',boxcar_kernel','same');
%   bconv = conv2(bconv',boxcar_kernel','same');
% 
%   filtered = gconv - bconv;
% else
%   filtered = gconv;
% end

% BPB 20100712: version of kernels listed as JWM above as slow but updates 
% in internal matlab make this version as fast or faster than previous.
% Tested on Matlab 2010a, ver 7.10.0.499

% BPB 20131008 working on three dimensional filtering
% JN 20130314 added normalization to convolutions to reduce edge problems

% gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');

    




if any(lobject)
   kx=-1.5*lobject(1):1:1.5*lobject(1);
    ky=-1.5*lobject(2):1:1.5*lobject(2);
   kz=-1.5*lobject(3):1:1.5*lobject(3);
   
   [KX,KY,KZ]=ndgrid(kx,ky,kz);
   KR=sqrt(KX.^2+KY.^2+KZ.^2);
   
      
    if lnoise == 0
        lnoise=.01;
    end
    
      
   noiseKernel=normalize(exp(-KR.^2./2/lnoise));

   smoothkernel=normalize(exp(-.5*(KX.^2/lobject(1)+KY.^2/lobject(2)+KZ.^2/lobject(3))));
kernal=smoothkernel-noiseKernel;
gconv = convnfft(image_array,kernal,'same');
    filtered = -gconv;
      
else
    
    if lnoise == 0
  gaussian_kernel = 1;
else      
  gaussian_kernel = normalize(...
    exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));
    end

    gaussian_kernel = repmat(gaussian_kernel,[length(gaussian_kernel),1,length(gaussian_kernel)]);
gaussian_kernel = gaussian_kernel.*permute(gaussian_kernel,[3,1,2]).*permute(gaussian_kernel,[2,3,1]);

ngconv=convnfft(ones(size(image_array)),gaussian_kernel,'same');
gconv = convnfft(image_array,gaussian_kernel,'same')./ngconv;

    filtered = gconv;
end

% Zero out the values on the edges to signal that they're not useful.     
lzero = max(lobject,ceil(5*lnoise));

% filtered(1:(round(lzero(1))),:,:) = 0;
% filtered((end - lzero(1) + 1):end,:,:) = 0;
% filtered(:,1:(round(lzero(2))),:) = 0;
% filtered(:,(end - lzero(2) + 1):end,:) = 0;
% filtered(:,:,1:(round(lzero(3)))) = 0;
% filtered(:,:,(end - lzero(3) + 1):end) = 0;
% JWM: I question the value of zeroing out negative pixels.  It's a
% nonlinear operation which could potentially mess up our expectations
% about statistics.  Is there data on 'Now centroid gets subpixel accuracy
% much more robustly'?  To choose which approach to take, uncomment one of
% the following two approaches.
% ERD: The negative values shift the peak if the center of the cntrd mask
% is not centered on the particle.

% % Approach A
% res = filtered;

% % Approach B
filtered(filtered < threshold) = 0;
res = filtered;