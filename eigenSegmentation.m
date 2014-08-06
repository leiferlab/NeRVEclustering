function Dim=eigenSegmentation(imageIn,hessFlag)

%based on Yi's segmentation algorithm, finds eigenvalues and eigen vectors
%of hessian matrix, then measures the local ordering of the eigen vectors
%as a way of producing image contrast.

if nargin==1
    hessFlag=1;
end


H=hessianMatrix(imageIn,4);
[Heig,HV]=hessianEig(H);
Hmax=max(Heig,[],3);
%Hmax=Hmax./max(Hmax(:));

se=strel('disk',5);
se=se.getnhood;
%se=true(n);
maxF=ordfilt2(Hmax,sum(se(:)),se);
Hmax=Hmax./maxF;
Hmax=bwmorph(Hmax>.3,'thin',Inf);
%Hmax=windowNormalize(Hmax,10);

% deltaEig=Heig(:,:,1)-Heig(:,:,2);
% X=HV{1,1}.*deltaEig;
% Y=HV{2,1}.*deltaEig;
% X(isnan(X))=0;
% Y(isnan(Y))=0;
% g=X.*smooth2a(X,3,3)+Y.*smooth2a(Y,3,3);
testData=[imageIn(:),reshape(Heig,[],ndims(imageIn))];
testData=zscore(testData);
testData(:,end)=3*testData(:,end);
[IDX,C,sumd,D] =kmeans(testData,2,'start','uniform','replicates',4);
D=bsxfun(@rdivide,D,sum(D,2));

DD=D(:,1)./(sum(D(:,1:2),2));

if mean(imageIn(DD>.9))<mean(imageIn(DD<.1))
     DD=1-DD;
end

IDX=DD>.5;
%sIDX=IDX>1.5;

if mean(imageIn(IDX>.9))<mean(imageIn(IDX<.1))
    IDX=~IDX;
end

Dim=reshape(IDX,size(imageIn));
if hessFlag
Dim=Dim & ~Hmax;
end

% DD=D(:,1)./(sum(D(:,1:2),2));
% Dim=reshape(DD,size(imageIn));
% if mean(imageIn(Dim>.9))<mean(imageIn(Dim<.1))
%     Dim=1-Dim;
% end


