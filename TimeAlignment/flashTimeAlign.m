function [idxOut,minABdist]=flashTimeAlign(flashA,flashB);
%function takes 2 series of flash times and attempts to match flash
%identities, flash A should be the more complete set, output will be the
%indicies of the first flashes that are determined to be aligned. 
AdistLin=pdist(flashA);
Adist=(squareform(AdistLin));
BdistLin=pdist(flashB);
Bdist=(squareform(BdistLin));

ABdist=pdist2(AdistLin',BdistLin');
[minABdist,minABdistIdx]=min(ABdist);
[minABdist,ib]=sort(minABdist,'ascend');
Asize=length(AdistLin);
Bsize=length(BdistLin);
for i=1:length(minABdistIdx)
    
idxA=indFromSquareIdx(minABdistIdx(ib(i)),Asize);
idxB=indFromSquareIdx(ib(i),Bsize);

idxOut(i,:)=[idxA(1),idxB(1)];
end

end

function idx=indFromSquareIdx(n,plength)
reDist=zeros(1,plength)';
reDist(n)=1;
reDist=triu(squareform(reDist));
[idx(1),idx(2)]=find(reDist);

end



