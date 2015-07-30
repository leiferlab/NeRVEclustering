function  [f, g]= gaussOverlap(A,B,scales,peaks,limit)


%Jeff N's version of the mex_GaussTransform from GMM mixturemodel code. I'm
%adding possibility to change scales and making it more transparent. The
%code takes point sets A and B, puts a gaussian at each position of with
%standard deviation given by scale, and then calculates the inner product.
%NOTE: this is not the full inner product as it skips the prefactor in
%front of the exponential term. This makes it so that the inner product of
%a gaussian with itself is 1. This is done to match the mex program
%exactly. Speed is now slightly faster than mex_GaussTransform

% inputs:   A pointset of scene, If A is n*3, scales and peaks are
% specified sepearately. if A is n*4, the 4th D is the scales, if A is n*5,
% 5th D is the peaks 
%           B pointset of model, same format as A;
%           scale is the variances of each point. if it is scalar, all
%           points will have the same varaince, if it is a cell, the first
%           cell will be the scales for the points in A, the second will be
%           the scales for the points in B

% outputs : f total overlap integrals
%            g gradient in energy of the points in A
   
Apoints=size(A,1);Bpoints=size(B,1);
nPoints=Apoints*Bpoints;


if nargin<5
    limit=inf;
end
if nargin<4
    if size(A,2)>4
        peaksA=A(:,5);
        peaksB=B(:,5);
        peaksA=peaksA/sum(peaksA);
        peaksB=peaksB/sum(peaksB);
        peaksMat=bsxfun(@times, peaksA,peaksB');

    else
    peaksMat=1/nPoints;
    end
else
            peaksA=peaks{1};peaksB=peaks{2};

    peaksMat=bsxfun(@times, peaksA,peaksB');
end
if ~iscell(scales)
        if size(A,2)>3
        scaleA=A(:,4)*scales;
        scaleB=B(:,4)*scales;
        varMat=bsxfun(@plus,scaleA.^2,(scaleB').^2);

        else
    varMat=2*scales.^2;
        end
    
else

     varMat=bsxfun(@plus,scales{1}.^2,(scales{2}').^2);
end

   

dmatx=bsxfun(@minus,A(:,1),B(:,1)');
dmaty=bsxfun(@minus,A(:,2),B(:,2)');
dmatz=bsxfun(@minus,A(:,3),B(:,3)');
dmat=sqrt(dmatx.^2+dmaty.^2+dmatz.^2);
%% overlap of gaussians is a gaussian of the distance with variance as the sum
% of variances

% gradient of that overlap (d expMat/ dA) is force

 %expMat=peaksMat./sqrt(2*pi*varMat).*exp(-1/2 * dmat.^2./varMat);
  expMat=peaksMat.*exp(-2 * dmat.^2./varMat);

    f=sum(sum(expMat));

expMat(dmat>limit)=0;
% fmatx=(1./(2*pi*varMat)).*expMat.*2.*dmatx;
% fmaty=(1./(2*pi*varMat)).*expMat.*2.*dmaty;
% fmatz=(1./(2*pi*varMat)).*expMat.*2.*dmatz;

fmatx=-(1./(varMat)).*expMat.*4.*dmatx;
fmaty=-(1./(varMat)).*expMat.*4.*dmaty;
fmatz=-(1./(varMat)).*expMat.*4.*dmatz;

g=[sum(fmatx,2),sum(fmaty,2),sum(fmatz,2)];





