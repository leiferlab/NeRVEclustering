function A=sparseTransitionCorr(A,W)

% sparseTransitionCorr approximates the correlation matrix of a sparse
% transition matrix consisting of only ones and zeros. It takes advantage
% of sparseness to quickly calculate mean and std and approximates the
% covariance without mean centering. Also accepts sparse weight matrix W.
% Weighted Version is NOT mean centered!!!

if ~issparse(A)
    A=sparse(A);
end

if nargin==1



Alength=size(A,2);
asum=sum(A);
amean=asum/Alength;
aSTD=sqrt((asum.*(1-amean).^2+(Alength-asum).*amean.^2)/Alength);
A=A'*A/Alength;
A=bsxfun(@rdivide,A,aSTD);
A=bsxfun(@rdivide,A,aSTD');

elseif nargin==2
    if ~issparse(W)
    W=sparse(W);
    end

aSTD=sqrt(sum(W,1));
A=(A.*W)'*A;
A=bsxfun(@rdivide,A,aSTD);
A=bsxfun(@rdivide,A,aSTD');

    
else
    error('Wrong number of inputs')
end

