function V=normalizeRange(X)
%The function takes a scalar or nd matrix X and normalize the matrix to a
%minimum of 0 and a maximum of 1. If all the values in the matrix are the
%same, they are all replaced with ones. 
if ~isempty(X)
L=numel(X);
XX=reshape(X,1,L);
if max(XX)==min(XX)
    V=ones(size(X));
else
V=(X-min(XX))/(max(XX)-min(XX));
end
else
    V=X;
end