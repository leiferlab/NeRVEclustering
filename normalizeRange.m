function V=NormalizeRange(X)
L=numel(X);
XX=reshape(X,1,L);
if max(XX)==min(XX)
    V=ones(size(X));
else
V=(X-min(XX))/(max(XX)-min(XX));
end