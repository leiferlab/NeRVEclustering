function [H_out,K_out,X_out,Y_out,Z_out]= MGCurveloop(varargin)
%same and MGcurve, but loops the rows of the XY and Z matrices
%JN


if nargin==3
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    nr=3;
else
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    nr=varargin{4};
end
    

if( all(X(1,:)==X(end,:)) && all(Y(1,:)==Y(end,:)) &&all(Z(1,:)==Z(end,:)))
    X(end,:)=[];
    Y(end,:)=[];
    Z(end,:)=[];
    
end

loop_size=size(X,1);
X=repmat(X,3,1);
Y=repmat(Y,3,1);
Z=repmat(Z,3,1);
Xs=smooth2b(X,2,2);
Ys=smooth2b(Y,2,2);
Zs=smooth2b(Z,2,2);


 [Xu,Xv]     =   gradient(Xs,nr);
[Xuu,Xuv]   =   gradient(Xu,nr);
[Xvu,Xvv]   =   gradient(Xv,nr);

[Yu,Yv]     =   gradient(Ys,nr);
[Yuu,Yuv]   =   gradient(Yu,nr);
[Yvu,Yvv]   =   gradient(Yv,nr);

[Zu,Zv]     =   gradient(Zs,nr);
[Zuu,Zuv]   =   gradient(Zu,nr);
[Zvu,Zvv]   =   gradient(Zv,nr);

% Reshape 2D ArraYs into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

Xu          =   [Xu Yu Zu];
Xv          =   [Xv Yv Zv];
Xuu         =   [Xuu Yuu Zuu];
Xuv         =   [Xuv Yuv Zuv];
Xvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(Xu,Xu,2);
F           =   dot(Xu,Xv,2);
G           =   dot(Xv,Xv,2);

m           =   cross(Xu,Xv,2);
p           =   sqrt(dot(m,m,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(Xuu,n,2);
M           =   dot(Xuv,n,2);
N           =   dot(Xvv,n,2);

[s,t] = size(Z);

% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
H = reshape(H,s,t);
% gm=smooth2a(gm,3,3);
% gc=smooth2a(gc,3,3);

X_out=X(loop_size:2*loop_size,:);
Y_out=Y(loop_size:2*loop_size,:);
Z_out=Z(loop_size:2*loop_size,:);
H_out=H(loop_size:2*loop_size,:);
K_out=K(loop_size:2*loop_size,:);

