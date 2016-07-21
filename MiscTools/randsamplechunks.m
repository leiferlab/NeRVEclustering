function output=randsamplechunks(n,p,chunkSize)

% works similar to rand perm, but instead of moving around all numbers, it
% permutes sets of numbers with size chunkSize
leftover=rem(p,chunkSize);
nChunks=floor(p/chunkSize);
y = randsample(1:n,nChunks+1,1)';
window=1:chunkSize;

allPoints=bsxfun(@plus, y(1:end-1),window)';
allPoints=allPoints(:);
%add leftovers
allPoints=[allPoints; (y(end)+(1:leftover))'];
output=mod(allPoints,n)+1;
