function h=colorplot(varargin)
x=varargin{1};
y=varargin{2};
    c=varargin{end};
x=x(:);
y=y(:);
c=c(:);

if nargin==3

h=surface([x,x]',[y,y]',[c,c]', ...
    'edgecolor','interp','facecolor','no','linew',2);
else
    z=varargin{3};
    z=z(:);
    
h=mesh([x,x]',[y,y]',[z,z]',[c,c]', ...
    'edgecolor','interp','facecolor','no','linew',2);    
    
end