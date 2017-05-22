function [idxOut,outputOffset]=flashTimeAlign2(flashA,flashB)
%function takes 2 sets of flash times and attempts to match flash
%identities, flash A should be equal to or larger than flashB, output will be the
%indicies of the first flashes that are determined to be the same, along
%with the indices of flashA. If length(flashB)==1, then the program will
%simply find the closes match in flashA to the time for flashB. 

%Inputs: flashA - a 1xN set of times corresponding to flashes on one camera
%        flashB - a 1xM set of times corresponding to flashes on another
%        camera. M should be <=N

%Outputs: indxOut - a 1x2 vector with indxOut(1) being the index of the
%flash in B that matches with indxOut(2) from flashB
%         outputOffset - the indices of flashA that are matched with an
%         index in flashB.

if length(flashB)==1
    [~,minDistIdx]=min(abs(flashA-flashB));
    idxOut=[1,minDistIdx];
    outputOffset=minDistIdx;
    
    return
end



best=Inf;

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



