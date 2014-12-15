function [imOut,Xw,Yw,Zw] = tpswarp3(img, outDim, moving, model)
%
% Description: Thin-Plane spline warping of the input image (img) to
% output image (imgw). The warping is defined by a set of reference points
% (Zp=[Xp Yp]) on the [img] image and a set of corresponding points (Zs)
% on the (imgw) image. The warping will translate Zp in img to Zs in imgw.
%
% Input:
% img - input image
% outDim - Output canvas ([W H])
% Zp - landmark in img
% Zs - landmark in canvas
% interp.method - interpolation mode('nearest', 'invdist', 'none')
% interp.radius - Max radius for nearest neighbor interpolation or
%                 Radius of inverse weighted interpolation
% interp.power - power for inverse weighted interpolation
%
% Output:
% imgw - warped image with no holes
% imgwr - warped image with holes
% map - Map of the canvas with 0 indicating holes and 1 indicating pixel
%
% Reference:
% F.L. Bookstein, "Principal Warps: Thin-Plate splines and the
% decomposition of deformations", IEEE Transaction on Pattern Analysis and
% Machine Intelligence, Vol. 11, No. 6, June 1989
%
% Author: Fitzgerald J Archibald
% Date: 07-Apr-09

% modified Dec 1 2014  for 3D data, Data must be 3D as the new kernal is
% hard coded. 
% modified Dec 12 2014 so image is not needed, You can input points instead
% and the points will move from the model space to moving space, this is
% done by interpolating Yw,Xw and Zw at the input points.

%% Initialization
NPs = size(moving,1); % number of landmark points
if isempty(outDim)
    if ismatrix(img)
        outDim=ceil(max([img; moving; model],[],1));
        outDim=outDim([2,1,3]);
    elseif ndims(img)==3;
        outDim=size(img);
    end
end

    
imgH = size(img,1); % height
imgW = size(img,2); % width
imgP=size(img,3);

% landmark in input
Xp = moving(:,2)';
Yp = moving(:,1)';
Zp=moving(:,3)';

% landmark in output (homologous)
Xs = model(:,2)';
Ys = model(:,1)';
Zs= model(:,3)';
%% Algebra of Thin-plate splines

% Compute thin-plate spline mapping [W|a1 ax ay] using landmarks
[wL]=conputeWl_3(moving);
wY = [Xs(:) Ys(:) Zs(:); zeros(4,3)]; % Y = ( V| 0 0 0)'   where V = [G] where G is landmark homologous (nx2) ; Y is col vector of length (n+3)
wW = (wL)\wY; %wW = inv(wL)*wY; (W|a1 ax ay)' = inv(L)*Y

% Thin-plate spline mapping (Map all points in the plane)
% f(x,y) = a1 + ax * x + ay * y + SUM(wi * U(|Pi-(x,y)|)) for i = 1 to n

[Yw, Xw,Zw]=tpsMap(wW, outDim,moving, NPs);


%% Warping

% input grid for warping
[X Y,Z] = ndgrid(1:outDim(1),1:outDim(2),1:outDim(3)); % HxW
Xw=reshape(Xw,outDim);
Yw=reshape(Yw,outDim);
Zw=reshape(Zw,outDim);
if ndims(img)==3
imOut=interp3(img, Yw,Xw,Zw, '*linear');
elseif ismatrix(img) && size(img,2)==3
    imOut=[interp3(Yw,img(:,1),img(:,2),img(:,3)) ...
        interp3(Xw,img(:,1),img(:,2),img(:,3)) ...
        interp3(Zw,img(:,1),img(:,2),img(:,3))];
        
    
else
    imOut=[];
end

% Nearest neighbor or inverse distance weighted based interpolation

return

%% [L] = [[K P];[P' 0]]
% np - number of landmark points
% (xp, yp) - coordinate of landmark points
function [wL]=computeWl(xp, yp, np)

rXp = repmat(xp(:),1,np); % 1xNp to NpxNp
rYp = repmat(yp(:),1,np); % 1xNp to NpxNp

wR = sqrt((rXp-rXp').^2 + (rYp-rYp').^2); % compute r(i,j)

wK = radialBasis(wR); % compute [K] with elements U(r)=r^2 * log (r^2)
wP = [ones(np,1) xp(:) yp(:)]; % [P] = [1 xp' yp'] where (xp',yp') are n landmark points (nx2)
wL = [wK wP;wP' zeros(3,3)]; % [L] = [[K P];[P' 0]]

return

function [wL]=conputeWl_3(x)

wR=squareform(pdist(x));
wK=wR;
wP = [ones(size(x,1),1) x]; % [P] = [1 xp' yp'] where (xp',yp') are n landmark points (nx2)
wL = [wK wP;wP' zeros(4)]; % [L] = [[K P];[P' 0]]


%% Mapping: f(x,y) = a1 + ax * x + ay * y + SUM(wi * U(|Pi-(x,y)|)) for i = 1 to n
% np - number of landmark points
% (xp, yp) - coordinate of landmark points
function [Xw, Yw,Zw]=tpsMap(wW, outDim, x,  np)

[X, Y ,Z] = ndgrid(1:outDim(1),1:outDim(2),1:outDim(3)); % HxW
X=X(:); % convert to 1D array by reading columnwise (NWs=H*W)
Y=Y(:); % convert to 1D array (NWs)
Z=Z(:);
NWs = length(X); % total number of points in the plane

% Mapping Algebra
wR=pdist2(x,[X Y Z] );
wK = (wR); % compute [K] with elements U(r)=r^2 * log (r^2)
wP = [ones(NWs,1) X(:) Y(:) Z(:)]'; % [P] = [1 x' y'] where (x',y') are n landmark points (nx2)
wL = [wK;wP]'; % [L] = [[K P];[P' 0]]

Xw  = wL*wW(:,1); % [Pw] = [L]*[W]
Yw  = wL*wW(:,2); % [Pw] = [L]*[W]
Zw = wL*wW(:,3);
return

%% k=(r^2) * log(r^2)
function [ko]=radialBasis(ri)

r1i = ri;
r1i(find(ri==0))=realmin; % Avoid log(0)=inf
ko = 2*(ri.^2).*log(r1i);

return
