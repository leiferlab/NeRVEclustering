function A=sparseTransitionCorr(A,W,m)

% sparseTransitionCorr approximates the correlation matrix of a sparse
% transition matrix consisting of only ones and zeros. It takes advantage
% of sparseness to quickly calculate mean and std and approximates the
% covariance without mean centering. Also accepts sparse weight matrix W.
% Weighted Version is NOT mean centered!!!

%inputs:
%   A: Sparse binary matrix of features, size MxN with N items and M
%   features
%   W: weight matrix with same size as A



if ~issparse(A)
    A=sparse(A);
end
if nargin==1
    W=[];
    m=1; 
end

%unweighted version
if isempty(W)
    Alength=size(A,2);
    asum=sum(A);
    amean=asum/Alength;
    aSTD=sqrt((asum.*(1-amean).^2+(Alength-asum).*amean.^2)/Alength);
    aSTD(asum==0)=1;
    %calculate correlation
    if nargin<3
        %calculate direction
        A=A'*A/Alength;
    else
        %calculate in chunks of 1000
        
        %build chunks
        win=1000;
        selectionRange=1:Alength;
        c=repmat(win,1,floor(Alength/win));
        c=[c mod(Alength,win)];
        selectionRange=mat2cell(selectionRange,1,c);
        Aout=cell(1,length(c));
        progressbar(0)
        %calculate for each chunk
        for iChunk=1:length(c)
            tic
            Atemp=A'*A(:,selectionRange{iChunk});
            %help sparsify a bit more by getting rid of elements with only
            %one match
            Atemp(Atemp<=m)=0;
            Aout{iChunk}=Atemp;
            %Aout=cat(2,Aout,Atemp);
            progressbar(iChunk/length(c));
        end
        A=cell2mat(Aout)/Alength;
        
    end
    A=bsxfun(@rdivide,A,aSTD);
    A=bsxfun(@rdivide,A,aSTD');
%weighted version, not mean subtracted, more computationally intensive
elseif nargin==2
    if ~issparse(W)
        W=sparse(W);
    end
    aSTD=sqrt(sum(W,1));
    A=(A.*W)'*A;
    A=bsxfun(@rdivide,A,aSTD);
    A=bsxfun(@rdivide,A,aSTD');
    
    
else
    error('Wrong number of inputs')
end

