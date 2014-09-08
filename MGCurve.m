function [H,K]= MGCurve(varargin)

%function takes a surface given by inputs X,Y,Z, and a smoothing kernal nr
%and calculates the mean and gaussian curvature everywhere on the surface.
%using the first and second fundamental forms. JN


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



 [Xu,Xv]     =   gradient(X,nr);
 Xu=smooth2a(Xu,2,2);
 Xv=smooth2a(Xv,2,2);
[Xuu,Xuv]   =   gradient(Xu,nr);
[Xvu,Xvv]   =   gradient(Xv,nr);

[Yu,Yv]     =   gradient(Y,nr);
 Yu=smooth2a(Yu,2,2);
 Yv=smooth2a(Yv,2,2);
[Yuu,Yuv]   =   gradient(Yu,nr);
[Yvu,Yvv]   =   gradient(Yv,nr);

[Zu,Zv]     =   gradient(Z,nr);
 Zu=smooth2a(Zu,2,2);
 Zv=smooth2a(Zv,2,2);
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

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);
