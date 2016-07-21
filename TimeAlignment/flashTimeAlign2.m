function [idxOut,outputOffset]=flashTimeAlign2(flashA,flashB);
%function takes 2 series of flash times and attempts to match flash
%identities, flash A should be the more complete set, output will be the
%indicies of the first flashes that are determined to be aligned. 
AdistLin=pdist(flashA);
Adist=(squareform(AdistLin));
BdistLin=pdist(flashB);
Bdist=(squareform(BdistLin));
best=Inf;
ABdist=pdist2(AdistLin',BdistLin');
[minABdist,minABdistIdx]=min(ABdist);
[minABdist,ib]=sort(minABdist,'ascend');
Asize=length(AdistLin);
Bsize=length(BdistLin);
for i=1:length(flashB)
    for j=1:length(flashA)
    offset=flashA(j)-flashB(i);
    ABdist=pdist2(flashA,flashB+offset);
    [minDistAll,minDistIdx]=min(ABdist);
    minDistAll=sum(minDistAll);
    if minDistAll<best
        best=minDistAll;
        idxOut=[i,j];
        outputOffset=minDistIdx;
        
    end
    
end
end


end

function idx=indFromSquareIdx(n,plength)
reDist=zeros(1,plength)';
reDist(n)=1;
reDist=triu(squareform(reDist));
[idx(1),idx(2)]=find(reDist);

end



