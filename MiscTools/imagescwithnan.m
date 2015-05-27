function [h hcb] = imagescwithnan(varargin)
% IMAGESC with NaNs assigning a specific color to NaNs
%imagescwithnan(c,cmap,nancolor,caxis);
%imagescwithnan(x,y,c,cmap,nancolor,caxis);

%# find minimum and maximum
    a=varargin{end-2};
    cm=varargin{end-1};
    nanclr=varargin{end};
        x=1:size(a,2);
    y=1:size(a,1);

if nargin>=5
x=varargin{1};
y=varargin{2};
end

amin=min(a(:));
amax=max(a(:));

%# size of colormap
n = size(cm,1);
%# color step
dmap=(amax-amin)/n;

%# standard imagesc
him = imagesc(x,y,a);
%# add nan color to colormap
colormap([nanclr; cm]);
%# changing color limits
if ~isnan(amin)||~isnan(amax)|| ~isnan(dmap)
caxis([amin-dmap amax]);

%# place a colorbar
hcb = colorbar;
%# change Y limit for colorbar to avoid showing NaN color
ylim(hcb,[amin amax])
end
if nargout > 0
    h = him;
end