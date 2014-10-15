function D = bwdist_jn(BW, scale)

%use the regular bwdist with a euclidean distance function with scaling
%factor in each dimension. By default, the scaling factor is [1,1,1];
if nargin==1
    D=bwdist(BW);
    return
end
oldsize=size(BW);
newsize=scale.*size(BW);

if ndims(BW)==3 && length(scale)==3
    %scale up, bwdist, scale down
    p = oldsize(1);
    q = oldsize(2);
    r = oldsize(3);
    
    p1 = newsize(1);
    q1 = newsize(2);
    r1 = newsize(3);
    
    [Y,X,Z] = meshgrid( (0:q-1)/(q-1), (0:p-1)/(p-1), (0:r-1)/(r-1)  );
    [YI,XI,ZI] = meshgrid( (0:q1-1)/(q1-1), (0:p1-1)/(p1-1), (0:r1-1)/(r1-1) );
    BW = interp3( Y,X,Z, (BW), YI,XI,ZI ,'*nearest');
    
D=bwdist(BW);
%D=image_resize(D,oldsize(1),oldsize(2),oldsize(3));

D=interp3(YI,XI,ZI,D,Y,X,Z,'*cubic');
elseif ismatrix(BW) && length(scale)==2
    BW=(imresize(BW,newsize));
    D=bwdist(BW);
D=imresize(D,oldsize);
else
    error('Binary Image must be 2D or 3D and scale must be length ndims(BW)')
end
