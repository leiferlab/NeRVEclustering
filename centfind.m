% PURPOSE:  Modification of Erik Dufresne's centroid.  Calculates the
% center of a bright spot using either the centroid method or by fitting to
% a 2D gaussian. Inspired by Grier & Crocker's feature for IDL, but greatly
% simplified and optimized for matlab
function [out,xout] = centfind(im,mx,sz,bpsz,opt)
% INPUT:
% im: image to process, particle should be bright spots on dark background
%     with little noise ofen an bandpass filtered brightfield image or a
%     nice fluorescent image
%
% mx: locations of local maxima to pixel-level accuracy from pkfnd.m
%
% sz: diamter of the window over which to average to calculate the centroid.
%     should be big enough
%     to capture the whole particle but not so big that it captures others.
%     if initial guess of center (from pkfnd) is far from the centroid, the
%     window will need to be larger than the particle size.  RECCOMMENDED
%     size is the long lengthscale used in bpass plus 2.
%
% bpsz: size input from bandpass must be an integer.  Set to zero if you
%       didn't use bandpass
%
% opt: structure containing optional options
%
%      interactive:  OPTIONAL INPUT set this variable to one and it will
%                    show you the image used to calculate each centroid,
%                    the pixel-level peak and the centroid
%
%      method:       OPTIONAL INPUT set this variable to zero to use
%                    gaussian fit method, set to one to use centroid method
%                    (Defaults to gauassian fit method).  If gaussian
%                    fitting fails for some reason will fall back to
%                    centroid
%
% NOTE:
%  - if pkfnd, and cntrd return more then one location per particle then
%    you should try to filter your input more carefully.  If you still get
%    more than one peak for particle, use the optional sz parameter in pkfnd
%  - If you want sub-pixel accuracy, you need to have a lot of pixels in
%    your window (sz>>1). To check for pixel bias, plot a histogram of the
%    fractional parts of the resulting locations
%  - It is HIGHLY recommended to run in interactive mode to adjust the
%    parameters before you analyze a bunch of images.
%
% OUTPUT:
% The output format depends on the 'mode' of calling.
%
% out: a N x 4 array containing the fit parameters
% ---- with mode "centroid"
%           out(:,1) is the x-coordinates of centroid
%           out(:,2) is the y-coordinates of centroid
% ---- with mode "Gaussian"
%           out(:,1) is the x-coordinates of Gaussian
%           out(:,2) is the y-coordinates of Gaussian
% ---- all modes
%           out(:,3) is the brightnesses from centroid calc
%           out(:,4) is the square radius of gyration from centroid calc
%
% xout: extra output for gaussian method, contains centroid data and gauss
%       parameters
% ---- with mode "centroid"
%           xout(:,[1,2,3,4,7]) = 0;
%           xout(:,5) is centroid x-coordinate
%           xout(:,6) is centroid y-coordinate
% ---- with mode "Gaussian"
%           xout(:,1) is the x width of Gaussian
%           xout(:,2) is the y width of Gaussian
%           xout(:,3) is the rotation angle of Gaussian
%           xout(:,4) is the offset of Gaussian
%           xout(:,5) is centroid x-coordinate
%           xout(:,6) is centroid y-coordinate
%           xout(:,7) is the peak intensity of the Gaussian
%
% Dependencies: Optimization Toolbox (for Gaussian Method Only)
%
% Author: Eric R. Dufresne, Yale University,
% Created: Feb 4 2005
%
%  5/2005 inputs diamter instead of radius
%  Modifications:
%  D.B. (6/05) Added code from imdist/dist to make this stand alone.
%  ERD (6/05) Increased frame of reject locations around edge to 1.5*sz
%  ERD 6/2005  By popular demand, 1. altered input to be formatted in x,y
%  space instead of row, column space  2. added forth column of output,
%  rg^2
%  ERD 8/05  Outputs had been shifted by [0.5,0.5] pixels.  No more!
%  ERD 8/24/05  Woops!  That last one was a red herring.  The real problem
%  is the "ringing" from the output of bpass.  I fixed bpass (see note),
%  and no longer need this kludge.  Also, made it quite nice if mx=[];
%  ERD 6/06  Added size and brightness output ot interactive mode.  Also
%   fixed bug in calculation of rg^2
%  JWM 6/07  Small corrections to documentation
%
% Version: 2.0
% Revisions: Weisshaar Group
% Version  Date       Author   Description
% 2.0      2009.09.30 CJI      Added option to use 2D gauss fit and changed
%                              some formating
% 2.1      2009.10.02 CJI      Fixed Gaussian Fitting with BPB's good
%                              guesses
% 2.2      2010.05.04 CJI      added bpsz stuff this needs to be fixed up
% 2.3      2012.07.31 BPB      added method of circular symmetry
%% Handle Inputs

if nargin==3
    interactive=0;
    method=0;
else
    if isfield(opt,'interactive')
        interactive = opt.interactive;
    else
        interactive = 0;
    end
    
    if isfield(opt,'method')
        method = opt.method;
    else
        method = 0;
    end
end

if sz/2 == floor(sz/2)
    warning('spt:centfind:sizeOdd','sz must be odd, like bpass');
    % BPB wonders, why does it need to be odd?
end

if isempty(mx)
    warning('spt:centfind:noInitialPts','there were no positions inputted into cntrd. check your pkfnd theshold')
    out=[];
    xout=[];
    return;
end

%% create mask - window around trial location over which to calculate the
r=(sz+1)/2;
m = 2*r;
x = 0:(m-1) ;
cent = (m-1)/2;
x2 = (x-cent).^2;
dst=zeros(m,m);
for i=1:m
    dst(i,:)=sqrt((i-1-cent)^2+x2);
end

msk=zeros([2*r,2*r]);
msk(dst < r)=1.0;

dst2=msk.*(dst.^2);
% ndst2=sum(sum(dst2));

[nr,nc]=size(im);

nozone = floor(sz/2)+bpsz;
%remove all potential locations within distance sz from edges of image
% ind=find(mx(:,2) > 1.5*sz & mx(:,2) < nr-1.5*sz);
mx=mx(mx(:,2) > nozone & mx(:,2) < nr-nozone,:);
% ind=find(mx(:,1) > 1.5*sz & mx(:,1) < nc-1.5*sz);
mx=mx(mx(:,1) > nozone & mx(:,1) < nc-nozone,:);

[nmx,~] = size(mx);

%inside of the window, assign an x and y coordinate for each pixel
xl=repmat(1:2*r,2*r,1);
yl=xl';

%% Center Finding
pts=[];
extra=[];
%loop through all of the candidate positions
for i=1:nmx
    %create a small working array around each candidate location, and apply the window function
    tmp=msk.*im((mx(i,2)-r+1:mx(i,2)+r),(mx(i,1)-r+1:mx(i,1)+r));
    
    
    %calculate the total brightness
    norm=sum(sum(tmp));
    
    %calculate the radius of gyration^2
    rg=(sum(sum(tmp.*dst2))/norm);
    
    switch method
        case {0,1} % 2D gaussian, centroid
            % centroid center
            xavg=sum(sum(tmp.*xl))./norm; %calculate the weigthed average x location
            yavg=sum(sum(tmp.*yl))./norm; %calculate the weighted average y location
        otherwise
    end
    
    
    % Fit to 2D gaussian or use centroid data to find center
    switch method
        case 0 % 2D gaussian of some type
            
            % using the centroid information, fit to a gauss2d
            b = [min(tmp,2),min(tmp,1)];
            a0(1) = mean(b(:));
            a0(2) = tmp(floor(yavg),floor(xavg));
            a0(3) = xavg;   % maybe use the center of the box 'r' here
            a0(4) = yavg;   % maybe use the center of the box 'r' here
            a0(5) = sqrt(rg);
            a0(6) = sqrt(rg);
            a0(7) = 0;
            
            % Do Fitting
            options = optimset('Display','none','Jacobian','off',...
                'MaxFunEvals',3000,'MaxIter',1000);
            [Param,~,~,exitflag,~]=lsqcurvefit(@gauss2D,a0,tmp,tmp,[],[],options);
            
            if exitflag < 1
                % if fitting fails output centroid method centers
                Param = zeros(7,1);
                xoffset = xavg;
                yoffset = yavg;
            else
                xoffset = Param(3);
                yoffset = Param(4);
            end
            
        case 1 % centroid method if method == 1
            Param = zeros(7,1);
            xoffset = xavg;
            yoffset = yavg;
            
            
        case 2 % circular symmetry 
            % BPB mentions, needs to be checked
            Param = zeros(7,1);
            [xavg, yavg] = radialcenter(tmp);
            xoffset = xavg;
            yoffset = yavg;
    end
    
    %concatenate it up
    pts = [pts,...
        [mx(i,1)+xoffset-r,...  %x center
        mx(i,2)+yoffset-r,...   %y center
        norm,...                %brightnesses
        rg]'];                  %radius of gyration (centroid)
    
    s2fwhm = 2*sqrt(2*log(2));
    % bpb asks: why save the centroid center if doing a gaussian fit?
    extra = [extra,...
        [s2fwhm*Param(5),...    %x width
        s2fwhm*Param(6),...     %y width
        Param(7),...            %rotation angle
        Param(1),...            %offset
        mx(i,1)+xavg-r,...      %x center (centroid)
        mx(i,2)+yavg-r,...      %y center (centroid)
        Param(2)]'];            %intensity
    
    %OPTIONAL plot things up if you're in interactive mode
    if interactive==1
        imagesc(tmp)
        axis image
        hold on;
        plot(xoffset,yoffset,'x')
        plot(xoffset,yoffset,'o')
        plot(r,r,'.')
        hold off
        title(['brightness ',num2str(norm),' size ',num2str(sqrt(rg))])
        pause
    end
    
    
end

%% Outputs
out=pts';
xout=extra';