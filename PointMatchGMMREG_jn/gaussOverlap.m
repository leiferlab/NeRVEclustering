function  [f, g]= gaussOverlap(A,B,scales,limit)


%Jeff N's version of the mex_GaussTransform from GMM mixturemodel code. I'm
%adding possibility to change scales and making it more transparent. The
%code takes point sets A and B, puts a gaussian at each position of with
%standard deviation given by scale, and then calculates the overlap. This
%is done very easily with Eqn 2 from http://www.tina-vision.net/docs/memos/2003-003.pdf


% inputs:   A pointset of scene
%           B pointset of model
%           scales standard devation of points, fornow, all the same.

% outputs : f total overlap integrals
%            g gradient in energy of the points in A
   
if nargin==3
    limit=inf;
end

dmat=pdist2(A,B);
dmatx=pdist2(A(:,1),B(:,1),@minus);
dmaty=pdist2(A(:,2),B(:,2),@minus);
dmatz=pdist2(A(:,3),B(:,3),@minus);
%% overlap of gaussians is a gaussian of the distance
expMat=-1/sqrt(4*pi*scales(1).^2)*exp(-1/4 * dmat.^2/scales(1).^2);

f=sum(sum(expMat));

expMat(dmat>limit)=0;
%% gradient of that overlap (d expMat/ dA) is force
fmatx=(1/(4*pi*scales(1).^2))*expMat.*2.*dmatx;
fmaty=(1/(4*pi*scales(1).^2))*expMat.*2.*dmaty;
fmatz=(1/(4*pi*scales(1).^2))*expMat.*2.*dmatz;

g=-[sum(fmatx,2),sum(fmaty,2),sum(fmatz,2)];





