function  r=sparseTransitionCorr(A)

% sparseTransitionCorr approximates the correlation matrix of a sparse
% transition matrix consisting of only ones and zeros. It takes advantage
% of sparseness to quickly calculate mean and std and approximates the
% covariance without mean centering


if ~issparse(A)
    A=sparse(A);
end


Alength=size(A,1);
asum=sum(A);
amean=asum/Alength;
aSTD=sqrt((asum.*(1-amean).^2+(Alength-asum).*amean.^2)/Alength);
aSTD=full(aSTD);
aCovariance=A*A'/Alength;
stdMat=(aSTD'*aSTD);
r=aCovariance./stdMat;
