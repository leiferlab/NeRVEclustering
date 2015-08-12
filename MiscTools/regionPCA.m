function [Dout,Vout,centroids]=regionPCA(cc)

% does pca on a region to approximate the size of objects. Input is
% connected component sturcutre for bwconncomp function.


if length(cc.ImageSize)==3        
    Dout=zeros(cc.NumObjects,3);
        Vout=zeros(cc.NumObjects,3,3);
centroids=zeros(cc.NumObjects,3);
    for i=1:length(cc.PixelIdxList);
        
        [x, y, z]=ind2sub(cc.ImageSize,cc.PixelIdxList{i});
        P=[x y z];
        P=bsxfun(@minus,P,mean(P));
        [v, d]=eig(P'*P);
        Dout(i,:)=d( logical(eye(3)));
        Vout(i,:,:)=v;
        centroids(i,:)=[mean(y),mean(x),mean(z)];
        
    end
else
        Dout=zeros(cc.NumObjects,2);
        Vout=zeros(cc.NumObjects,2,2);
centroids=zeros(cc.NumbObjects,2);
        for i=1:length(cc.PixelIdxList);
        

        [x, y]=ind2sub(cc.ImageSize,cc.PixelIdxList{i});
        P=[x y];
        P=bsxfun(@minus,P,mean(P));
        [v, d]=eig(P'*P);
        Dout(i,:)=d( logical(eye(2)));
        Vout(i,:,:)=v;
                centroids(i,:)=[mean(y),mean(x)];

        
    end
end

