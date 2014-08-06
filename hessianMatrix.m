function [H,D]=hessianMatrix(A,h)
% function computes the hessian matrix of A at all points answer is output
% as a cell array. Also computes the divergence,D. input matrix A can be
% either 2 or 3D, h is the smoothing distance over which to take descrete
% differences. The outputs are in the form of cell arrays.

if nargin==1
    h=1;
end
if ndims(A)==3
[Ay,Ax,Az]=gradient(A,h);
[Axy,Axx,Axz]=gradient(Ax,h);
[Ayy,Ayx,Ayz]=gradient(Ay,h);
[Azy,Azx,Azz]=gradient(Az,h);

H{1,1}=Axx; H{1,2}=Ayx; H{1,3}=Azx;
H{2,1}=Axy; H{2,2}=Ayy; H{2,3}=Azy;
H{3,1}=Axz; H{3,2}=Ayz; H{3,3}=Azz;

D{1,1}=Ax; D{1,2}=Ay; D{1,3}=Az;

elseif ndims(A)==2
[Ay,Ax]=gradient(A,h);
[Axy,Axx]=gradient(Ax,h);
[Ayy,Ayx]=gradient(Ay,h);


H{1,1}=Axx; H{1,2}=Ayx; 
H{2,1}=Axy; H{2,2}=Ayy;
    
D{1,1}=Ax; D{1,2}=Ay;
end


