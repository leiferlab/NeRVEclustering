function tform =makeAffine3d(x,y,w)

% function makes tform structure for 3d set of points find the affine
% transform that takes the points x to the points y. optional input w is
% the weighting for each of the coordinates
x=[x ones(size(x(:,1)))];
A1=reshape([x(:) zeros(size(x(:))) zeros(size(x(:)))]', [],4);
A=[A1 circshift(A1,[1,0]) circshift(A1,[2,0])];
if nargin==2
Ainv=A\y(:);
else
    
end

y=[y ones(size(y(:,1)))];
Ainv=reshape(Ainv,4,3);
Ainv=[Ainv ,[0;0;0;1]];
Ainv=(y'/x')';
Ainv(:,end)=[0 0 0 1];
tform=affine3d(Ainv);
    

function [G,F] = wpinv(A,W)
  
% WPINV - Compute weighted pseudoinverse.
% 
%  [G,F] = wpinv(A,W)
% 
% Computes the optimal solution x = Gy + Fx0 to the least squares
% problem
% 
%  min ||W(x-x0)||  subj. to  Ax = y
%   x
%
% See also PINV.
  
  [m,n] = size(A);
  if nargin < 2
    W = eye(n);
  end
  
  % Thesis, Lemma B.1
  G = inv(W)*pinv(A*inv(W));
  F = eye(n) - G*A;