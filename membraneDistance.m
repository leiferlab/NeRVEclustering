function mDist=membraneDistance(P1,P2,V1)
%P1, P2 in 'patch' format
% V1 is face normals
% Distance from 2 to 1 (sign convention might have issues)
if nargin==2
    if ndims(P1)==3
    if  ~isfield(P1,'vertices')
        V1=vertnorm_dw(P1.faces,P1.vertices);
    else
        [~,~,V1]=surfArea(P1(:,:,1),P1(:,:,2),P1(:,:,3)); 
    end
    
    else
        error('Not enough Inputs')
    end
    
end


V1=reshape(V1,[],3);
P1=reshape(P1,[],3);
P2=reshape(P2,[],3);
[dist,I] = pdist2(P1,P2,'euclidean','Smallest',1);
V2=V1(I,:);
dM=P2-P1(I,:);

mDist=dot(dM,-V2,2);


