function feature_out=oneHotNeuron(features,n_outcomes)
% oneHotNeuron turns an positive integer matrix into a binary one-hot-encoded
% matrix.
%Inputs: 
%   features:
%          an n x k matrix of integers. 0's will be ignored
%   n_outcomes:
%          a 1 x k vector corresponding with each element as the maximum
%          possible number for each element of feature
%Outputs:
%   feature_out:
%           a n x sum(k) vector of zeros and ones, will be sparse
index_add=[0,cumsum(n_outcomes)];

if any(features(:))
    validPoints=features(:)>0;
    features=bsxfun(@plus,features ,index_add(1:end-1));
    
    %add indexing offset for sample
    startPos=(1:size(features,1));
    startPos=repmat(startPos,size(features,2),1)';
    transitionIdxX=startPos(validPoints);
    transitionIdxY=features(validPoints);
    
    
    %build binary sparce matrix, referred to as a transition matrix
    nTransitions=length(transitionIdxY);
    feature_out=sparse(transitionIdxX,transitionIdxY,ones(1,nTransitions),...
        size(features,1),max(index_add),nTransitions);
else
    if isempty(features)
        feature_out=[];
    else
        feature_out=sparse(size(features,1),max(index_add));
    end
end




