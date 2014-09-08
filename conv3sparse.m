function [convOut,offset]= conv3sparse(varargin)
% conv3sparse takes the convolution of input{1}  and {2} using a rolling
% sum. It excludes zeros in the original. This was in hopes of speeding
% things up, but convnfft is faster so I've been using that instead.

%JN



    A=varargin{1};
    B=varargin{2};
if nargin==2
    shape='same';
    offset=0;
elseif nargin==3
    offset=0;
    shape=varargin{3};
else
    offset=varargin{4};
    shape=varargin{3};
end


    


if length(size(B))<4
n_pts=find(abs(A)>(max(abs(A(:)))/200));
n_val=A(n_pts);

[nx,ny,nz]=ind2sub(size(A),n_pts);
convOut=padarray(zeros(size(A)),(size(B)-size(B)./size(B)),'post');

for i=1:length(n_val)

    rangeX=nx(i):nx(i)+size(B,1)-1;
    rangeY=ny(i):ny(i)+size(B,2)-1;
    rangeZ=nz(i):nz(i)+size(B,3)-1;

    convOut(rangeX,rangeY,rangeZ)=convOut(rangeX,rangeY,rangeZ)+n_val(i)*B;

end
else
    n_pts=find(A);
n_val=A(n_pts);

[nx,ny,nz]=ind2sub(size(A),n_pts);
convOut=padarray(zeros(size(A)),(size(B)-size(B)./size(B)),'post');

for i=1:length(n_val)
    rangeX=nx(i):nx(i)+size(B,1)-1;
    rangeY=ny(i):ny(i)+size(B,2)-1;
    rangeZ=nz(i):nz(i)+size(B,3)-1;

    psfSlice=(nz(i)-offset);
    
    if psfSlice>size(B,4)
        psfSlice=size(B,4);
    elseif psfSlice<1
        psfSlice=1;
    end
    
    convOut(rangeX,rangeY,rangeZ)=convOut(rangeX,rangeY,rangeZ)+n_val(i)*B(:,:,:,psfSlice);

end
end

if strcmp(shape,'full')
elseif strcmp(shape,'same')
    trim(1)=(size(B,1)-1)/2;
    trim(2)=(size(B,2)-1)/2;
    trim(3)=(size(B,3)-1)/2;
    x=trim(1)+(1:size(A,1));
    y=trim(2)+(1:size(A,2));
    z=trim(3)+(1:size(A,3));
    
    if all(~mod(trim,1))
        convOut=convOut(x,y,z);
    else
        [X,Y,Z]=meshgrid(x,y,z);
    convOut=interp3(convOut,X,Y,Z,'*spline');
    end
    
end


