function [ Heig ,HeigVec ] = hessianEig( H ,Hmask)

%UNTITLED Summary of this function goes here
%   Calculates eigenvalues of the hessian matrix H

if nargin==1
    Hmask=ones(size(H{1,1}));
end
Hsize=size(H{1,1});
Hind=find(Hmask);
if numel(Hsize)==3
    [Hmaskx,Hmasky,Hmaskz]=ind2sub(Hsize,Hind);
    
    Hlin=cellfun(@(x) reshape(x,1,1,[]),H,'UniformOutput',false);
    Hlin=cell2mat(Hlin);
    Heig=zeros(Hsize(1),Hsize(2),Hsize(3),3);
    for iPix=1:length(Hmaskx)
        [V,D]=eig(Hlin(:,:,Hind(iPix)));
        
        
        Heig(Hmaskx(iPix),Hmasky(iPix),Hmaskz(iPix),:)=...
            D(eye(3)>0);
        
    end  
elseif numel(Hsize)==2
        [Hmaskx,Hmasky]=ind2sub(Hsize,Hind);
    
    Hlin=cellfun(@(x) reshape(x,1,1,[]),H,'UniformOutput',false);
    Hlin=cell2mat(Hlin);
    Heig=zeros(Hsize(1),Hsize(2),2);
    for iPix=1:length(Hmaskx)
            [V,D]=eig(Hlin(:,:,Hind(iPix)));
           Heig(Hmaskx(iPix),Hmasky(iPix),:)=D(eye(2)>0);

    end  
end

if nargout==2
    
    
end


end
