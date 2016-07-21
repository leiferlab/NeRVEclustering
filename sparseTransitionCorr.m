function A=sparseTransitionCorr(A,W,m)

% sparseTransitionCorr approximates the correlation matrix of a sparse
% transition matrix consisting of only ones and zeros. It takes advantage
% of sparseness to quickly calculate mean and std and approximates the
% covariance without mean centering. Also accepts sparse weight matrix W.
% Weighted Version is NOT mean centered!!!

if ~issparse(A)
    A=sparse(A);
end
if nargin==1
    W=[];
end

if isempty(W)



Alength=size(A,2);
asum=sum(A);
amean=asum/Alength;
aSTD=sqrt((asum.*(1-amean).^2+(Alength-asum).*amean.^2)/Alength);
aSTD(asum==0)=1;
if nargin<3
A=A'*A/Alength;
else
    win=1000;
    selectionRange=1:Alength;
    c=repmat(win,1,floor(Alength/win));
    c=[c mod(Alength,win)];
    selectionRange=mat2cell(selectionRange,1,c);
    %Aout=spalloc(Alength,Alength,round((Alength/70)^2));
    %Aout=[];
    Aout=cell(1,length(c));
   progressbar(0)
    for iChunk=1:length(c)
        tic
        Atemp=A'*A(:,selectionRange{iChunk});
        Atemp(Atemp==m)=0;
        Aout{iChunk}=Atemp;
        %Aout=cat(2,Aout,Atemp);
        progressbar(iChunk/length(c));
    end
    A=cell2mat(Aout)/Alength;
    %A=Aout;
    
end
A=bsxfun(@rdivide,A,aSTD);
A=bsxfun(@rdivide,A,aSTD');

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

