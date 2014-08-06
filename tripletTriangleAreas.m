function A=tripletTriangleAreas(Apos)

[npts,ndims]=size(Apos);

%create distance matrix
Alldis=pdist(Apos);
Alldis=squareform(Alldis);

%all possible combinations of idx
[Yidx,Xidx,Zidx]=meshgrid(1:npts,1:npts,1:npts);

%look up edge lengths
edge1=Alldis(sub2ind(size(Alldis),Xidx(:),Yidx(:)));
edge2=Alldis(sub2ind(size(Alldis),Yidx(:),Zidx(:)));
edge3=Alldis(sub2ind(size(Alldis),Zidx(:),Xidx(:)));

%perimeter
S=edge1+edge2+edge3;

%Area via semiperimter
A=sqrt(S.*(S-edge1).*(S-edge2).*(S-edge3));
A=reshape(A,size(Yidx));

