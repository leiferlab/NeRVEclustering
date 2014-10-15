function imOut=watershedFilter(imIn,hmin,conn)
% simple standard watershed filtering

if nargin<3
        if ndims(imIn)==2
    conn=4;
    elseif ndims(imIn)==3
        conn=6;
        end
end

if nargin==1;
    hmin=0;

end

dIm=-bwdist(~imIn);
dIm(dIm==0)=Inf;
dIm=imhmin(dIm,hmin);
waterIm=watershed(dIm,conn);
imIn(~waterIm)=0;
imOut=imIn;