function [ Heig ,HeigVec ] = hessianEig( H ,Hmask)

%UNTITLED Summary of this function goes here
%   Calculates eigenvalues of the hessian matrix H

if nargin==1
    Hmask=true(size(H{1,1}));
else
    Hmask=Hmask>0;
end
Hsize=size(H{1,1});
Hind=find(Hmask);
if numel(Hsize)==3
  %  [Hmaskx,Hmasky,Hmaskz]=ind2sub(Hsize,Hind);
    
    Hlin=cellfun(@(x) reshape(x,1,1,[]),H,'UniformOutput',false);
    Hlin=cell2mat(Hlin);
    Heig=zeros(Hsize(1),Hsize(2),Hsize(3),3);
        [V,D]=eig3D(Hlin(:,:,Hind));
        
        
%            Heig(:,:,:,1)=reshape(D(1,1,:),Hsize);
%            Heig(:,:,:,2)=reshape(D(1,2,:),Hsize);
%            Heig(:,:,:,3)=reshape(D(1,3,:),Hsize);
for i=1:3
    htemp=zeros(Hsize);
    htemp(Hmask)=D(1,i,:);
    Heig(:,:,:,i)=htemp;
    if nargout==2
    for j=1:3
        
           HeigVec{i,j}=reshape(V(i,j,:),Hsize);
    end
    end
    
end




        
      
elseif numel(Hsize)==2
        [Hmaskx,Hmasky]=ind2sub(Hsize,Hind);
    
    Hlin=cellfun(@(x) reshape(x,1,1,[]),H,'UniformOutput',false);
    Hlin=cell2mat(Hlin);
    Heig=zeros(Hsize(1),Hsize(2),2);
    
    
            [V,D]=eig2D(Hlin(:,:,Hind));
         
            
for i=1:2
    htemp=zeros(Hsize);
    htemp(Hmask)=D(1,i,:);
    Heig(:,:,i)=htemp;
    if nargout==2
    for j=1:2
           htemp2=zeros(Hsize);
           htemp2(Hind)=V(i,j,:);
           
           HeigVec{i,j}=htemp2;
    end
    end
    
end
            
end

if nargout==2
    
    
end


end

function [vec,v]=eig2D(A)
%finds eigenvectors and eigen values of large set of 2x2 matricies, in put
%A is a 2x2xn matrix, output is v, a 1x2xn matrix of eigenvalues in
%descending order, vec is a 2x2xn matrix with corresponding eigenvectors.
a=A(1,1,:);b=A(1,2,:);c=A(2,1,:);d=A(2,2,:);
vec=zeros(2,2,size(A,3));
T=a+d;D=a.*d-b.*c;
desc=sqrt(T.^2-4*(D));
desc=real(desc);
v=[(T+desc)/2,(T-desc)/2];
vec(:,1,:)=bsxfun(@rdivide,[ones(size(b)),-b./(a-v(1,1,:))],(1+(b./(a-v(1,1,:))).^2).^.5);
vec(:,2,:)=bsxfun(@rdivide,[ones(size(b)),-b./(a-v(1,2,:))],(1+(b./(a-v(1,2,:))).^2).^.5);
end

function [vec,v]=eig3D(A)
% shamelessly copied from wikipedia

% Given a real symmetric 3x3 matrix A, compute the eigenvalues
 
p1 = A(1,2,:).^2 + A(1,3,:).^2 + A(2,3,:).^2;

   q = (A(1,1,:)+A(2,2,:)+A(3,3,:))/3;
   p2 = (A(1,1,:) - q).^2 + (A(2,2,:) - q).^2 + (A(3,3,:) - q).^2 + 2 * p1;
   p = sqrt(p2 / 6);
   B = bsxfun(@times,(1 ./ p),  (A - bsxfun(@times,q,eye(3)))) ;     
   r = B(1,1,:).*B(2,2,:).*B(3,3,:) + 2*B(1,2,:).*B(1,3,:).*B(2,3,:)...
       -B(1,1,:).*B(2,3,:).^2-B(2,2,:).*B(1,3,:).^2-B(3,3,:).*B(1,2,:).^2; %det of B
 r=r/2;
   % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
   % but computation error can leave it slightly outside this range.
   if (r <= -1) 
      phi = pi / 3;
   elseif (r >= 1)
      phi = 0;
   else
      phi = acos(r) / 3;
   end
 
   % the eigenvalues satisfy eig3 <= eig2 <= eig1
   eig1 = q + 2 * p .* cos(phi);
   eig3 = q + 2 * p .* cos(phi + (2*pi/3));
   eig2 = 3 * q - eig1 - eig3 ;    % since trace(A) = eig1 + eig2 + eig3

v=[eig1,eig2,eig3];
vec=zeros(size(A));
   for iEig=1:3
   % find eigenvectors using cramers rule to solve systems after setting z
   % component to 1, then normalize.
   D=(A(1,1,:)-eig1).*(A(2,2,:)-eig1)-A(1,2,:).^2;
   DX=(A(1,1,:)-eig1).*(-A(2,3,:))-(A(1,2,:).*(-A(1,3,:)));
   DY=(A(2,2,:)-eig1).*(-A(1,3,:))-(A(1,2,:).*(-A(2,3,:)));
   
   
   x1=DX./D;
   y1=DY./D;
   z1=ones(size(x1));
   
   L=sqrt(x1.^2+y1.^2+z1.^2);
   x1=x1./L;
   y1=y1./L;
   z1=z1./L;
   
   vec(:,iEig,:)=[x1;y1;z1];
   end
   

   
end
