function [P,k]=powerSpectrum(f,x,normlaizeFlag)
%function takes input f and stepsize x to calculate the corresponding power
%specturm P for each row and the kspace index k
if nargin==1;
    x=1;
    normlaizeFlag=0;
elseif nargin==2
        normlaizeFlag=0;

end

N=size(f,2);
F=fft(f,[],2);
%F=fftshift(F,2);
P=F.*conj(F);
P=P(:,1:N/2+1);
P(:,2:end)=2*P(:,2:end);
P=P*x/N;
if  normlaizeFlag==1;
P=bsxfun(@rdivide,P,sum(P,2))*(x*N);
end
%k=1./x;
k=0:1/(x*N):1/(2*x);
