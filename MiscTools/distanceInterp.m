function xyzout=distanceInterp(xyzs,n)
%interpolates curve to create n evenly spaced points along the curve
%descirbed by points xyzs
% if one input, just reinterp points to evenly space
if nargin==1
    n=size(xyzs,1);
end

d=size(xyzs,2);
            xyzout=zeros(n,d);
            ds=sqrt(sum(diff(xyzs).^2,2));
s=[0;cumsum(ds)];
Lserach=linspace(0,max(s),n);
for iDim=1:d;
    xyzout(:,iDim)=interp1(s,xyzs(:,iDim),Lserach,'spline');
end
