function output=randpermchunks(n,chunkSize)

% works similar to rand perm, but instead of moving around all numbers, it
% permutes sets of numbers with size chunkSize
chunkStart=ceil(rand*chunkSize);


chunks=bsxfun(@plus,1:chunkSize,(0:chunkSize:n-2*chunkSize-1)')+chunkStart;
chunk1={1:chunkStart};
chunkEnd={(max(chunks(:))+1):n};
chunks=mat2cell(chunks,ones(1,size(chunks,1)),size(chunks,2));
chunks=[chunk1;chunks;chunkEnd];
chunkOrder=randperm(length(chunks));

chunks=chunks(chunkOrder);
output=cell2mat(chunks');