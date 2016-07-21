function d=L2Dist(X,Y)
%trying to make a faster pdist2 hardcoded for euclidean L2 distance
% X=permute(X, [1 3 2]);
% Y=permute(Y,[3 1 2]);
% d=bsxfun(@minus,X,Y);
% d=sqrt(sum(d.^2,3));

        d = pdist2mex(X',Y','euc',[],[],[]);

