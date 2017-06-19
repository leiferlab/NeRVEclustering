function [Tv,Nv,Bv]=WormCurveCoordinates(centerline)
% make coordinate system around the worm by building the tangent, normal,
% and binormal vectors around a 3D curve. 
% Input:
%   centerline, an N*3 matrix of x,y,z coordinates
cl_size=size(centerline);
if numel(cl_size)==2
    cl_size(3)=1;
end
  Tv=zeros(cl_size(1),3,cl_size(3));
    Bv=Tv; Nv=Tv;
    % for each centerline, make the Tangent, Normal, and Binormal vectors
    for iSlice=1:cl_size(3)
        % T, N and B are first made in 2D,
        current_cl=centerline(:,:,iSlice);
        T=normr(gradient(current_cl',5)');
        N=[T(:,2) -T(:,1)];
        B=T(:,1).*N(:,2)-T(:,2).*N(:,1);
        N=bsxfun(@times, N,sign(B));
        B=sign(B);
        % and then the Zcomponent is added
        T=[T zeros(size(current_cl(:,1)))];
        N=[N zeros(size(current_cl(:,1)))];
        B=[zeros(size(current_cl(:,1))) zeros(size(current_cl(:,1))) B];
        %compile results into vectors containing TBN for each centerline
        Tv(:,:,iSlice)=T;
        Nv(:,:,iSlice)=N;
        Bv(:,:,iSlice)=B;
        
    end
    %fix any sign flipping in the binormal and normal vector. Make the
    %binormal always be positive, tangent points from tail to head, normal
    %fits in according the right hand rule
    signVector=sign(Bv(:,3,:));
    Bv=bsxfun(@times,Bv,signVector);
    Nv=bsxfun(@times,Nv,signVector);