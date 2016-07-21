function convOut = convnChunkfft(A,B,zchunk)

% calls convnfft in chunks and adds them up. This was in attempt to reduce
% the maximum memory usage, but we've since transitioned to using cetus and
% 64bit machines.


zchunk=max(round((450000-4*numel(B))/3/numel(A(:,:,1))),1);
zchunk=5;
if zchunk>size(A,3)
    zchunk=size(A,3);
end
A(isnan(A))=0;
convOut=padarray(zeros(size(A)),(size(B)-size(B)./size(B)),'post');
zsubs=1:zchunk:size(A,3);
if zsubs(end)~=size(A,3);
    zsubs=[zsubs,size(A,3)];
end


for i= 1:length(zsubs)-1
C = convnfft(A(:,:,zsubs(i):zsubs(i+1)), B, 'full');
convOut(:,:,zsubs(i):zsubs(i+1)+size(B,3)-1)=convOut(:,:,zsubs(i):zsubs(i+1)+size(B,3)-1)+C;
end

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
    